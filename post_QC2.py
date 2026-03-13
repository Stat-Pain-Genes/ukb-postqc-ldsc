#!/usr/bin/env python3
"""
post_QC.py — Post-filter REGENIE GWAS output
Filtres: EMAC, MAF, P-value, HWE (via fichier externe PLINK)

⚠️  IMPORTANT — Filtre HWE :
    Le test HWE ne peut PAS être calculé correctement depuis les summary stats
    REGENIE (comptes de génotypes observés absents). Deux options :

    Option A — Fichier HWE pré-calculé par PLINK (recommandé) :
        plink  --bfile mydata --hardy --out hwe_stats   # PLINK1 → .hwe
        plink2 --pfile mydata --hardy --out hwe_stats   # PLINK2 → .hardy
        python3 post_QC.py --in gwas.regenie.gz --hwe-file hwe_stats.hwe --hwe-minp 1e-12

    Option B — Désactiver HWE (si déjà filtré en amont) :
        python3 post_QC.py --in gwas.regenie.gz --no-hwe

Format REGENIE attendu (espace ou tab séparé):
  CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA

Usage:
  python3 post_QC.py --in gwas.regenie.gz --out gwas_filtered.tsv.gz --no-hwe
  python3 post_QC.py --in gwas.regenie.gz --out gwas_filtered.tsv.gz \
                     --hwe-file hwe_stats.hwe --hwe-minp 1e-12 \
                     --emac-min 100 --maf-min 0.001
"""
import sys, gzip, argparse, math
from typing import Optional, Dict

# ─────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────
def open_any(path: str):
    if path == "-" or path is None:
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def open_out(path: str):
    if path == "-":
        return sys.stdout
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")

def parse_float_safe(x: str) -> Optional[float]:
    try:
        v = float(x)
        return None if math.isnan(v) or math.isinf(v) else v
    except Exception:
        return None

# ─────────────────────────────────────────────────────────────
# HWE exact test (Wigginton et al., 2005)
# Utilisé uniquement quand les vrais comptes de génotypes sont fournis.
# ─────────────────────────────────────────────────────────────
def hwe_pvalue(obs_hom1: int, obs_het: int, obs_hom2: int) -> float:
    obs_homc = min(obs_hom1, obs_hom2)
    obs_homo = max(obs_hom1, obs_hom2)
    n        = obs_homc + obs_homo + obs_het
    if n == 0:
        return 1.0
    rare = 2 * obs_homc + obs_het
    mid  = int(rare * (2 * n - rare) / (2 * n))
    if (rare % 2) ^ (mid % 2):
        mid += 1
    prob_mid = 1.0
    tail     = prob_mid
    probL = prob_mid
    j = mid
    while j > 0:
        probL *= j * (j - 1) / (4.0 * (((rare - j) // 2) + 1) * ((((2*n - rare) - j) // 2) + 1))
        if probL <= prob_mid + 1e-15:
            tail += probL
        j -= 2
    probR = prob_mid
    k = mid
    while k <= rare - 2:
        probR *= (((rare - k) // 2) * (((2*n - rare) - k) // 2) * 4.0) / ((k + 2) * (k + 1))
        if probR <= prob_mid + 1e-15:
            tail += probR
        k += 2
    return min(1.0, tail)

# ─────────────────────────────────────────────────────────────
# Chargement fichier HWE externe (PLINK1 .hwe ou PLINK2 .hardy)
# ─────────────────────────────────────────────────────────────
def load_hwe_file(path: str) -> Dict[str, float]:
    """
    Charge un fichier HWE PLINK1 (.hwe) ou PLINK2 (.hardy).
    Retourne un dict  variant_id -> hwe_pvalue.

    Format PLINK1 .hwe (espace-separe, filtre sur TEST == "ALL") :
        CHR  SNP  TEST  A1  A2  GENO        O(HET)  E(HET)  P
        1    rs1  ALL   A   G   100/50/200  0.31    0.33    0.45

    Format PLINK2 .hardy (tab-separe) :
        #CHROM  ID  REF  ALT  HOM_REF_CT  HET_CT  HOM_ALT_CT  TWO_HIT_CT  O(HET)  E(HET)  P
    """
    hwe_dict: Dict[str, float] = {}
    n_loaded = 0
    n_skipped = 0

    try:
        fh = open_any(path)
    except FileNotFoundError:
        sys.stderr.write(f"Fichier HWE introuvable : {path}\n")
        sys.exit(1)

    header_line = fh.readline().rstrip("\n")

    # Detection du format
    # PLINK2 : commence par '#CHROM'
    # PLINK1 : commence par 'CHR' (sans #)
    is_plink2 = header_line.lstrip().startswith("#")

    if is_plink2:
        # PLINK2 .hardy — tab-separe
        sep_hwe = "\t"
        cols = [c.lstrip("#").strip() for c in header_line.split(sep_hwe)]
        id_col = next((i for i, c in enumerate(cols) if c in ("ID",)), None)
        p_col  = next((i for i, c in enumerate(cols)
                       if c in ("P", "P_MIDP", "P_EXACT")), None)

        if id_col is None or p_col is None:
            sys.stderr.write(f"Format PLINK2 .hardy non reconnu. Colonnes: {cols}\n")
            sys.exit(1)

        sys.stderr.write(f"  Format detecte : PLINK2 .hardy\n")
        sys.stderr.write(f"  Colonnes : ID={cols[id_col]}  P={cols[p_col]}\n")

        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split(sep_hwe)
            if len(parts) <= max(id_col, p_col):
                n_skipped += 1
                continue
            snp_id = parts[id_col].strip()
            p_val  = parse_float_safe(parts[p_col])
            if snp_id and p_val is not None:
                if snp_id not in hwe_dict or p_val < hwe_dict[snp_id]:
                    hwe_dict[snp_id] = p_val
                n_loaded += 1
            else:
                n_skipped += 1

    else:
        # PLINK1 .hwe — espace-separe
        cols = header_line.split()
        try:
            snp_i  = cols.index("SNP")
            test_i = cols.index("TEST")
            p_i    = cols.index("P")
        except ValueError:
            sys.stderr.write(f"Format PLINK1 .hwe non reconnu. Colonnes: {cols}\n")
            sys.exit(1)

        sys.stderr.write(f"  Format detecte : PLINK1 .hwe\n")
        sys.stderr.write(f"  Colonnes : SNP={cols[snp_i]}  TEST={cols[test_i]}  P={cols[p_i]}\n")

        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split()
            if len(parts) <= max(snp_i, test_i, p_i):
                n_skipped += 1
                continue
            # PLINK1 genere une ligne par type de test (ALL, AFF, UNAFF)
            # On garde uniquement TEST == "ALL"
            if parts[test_i].strip() != "ALL":
                continue
            snp_id = parts[snp_i].strip()
            p_val  = parse_float_safe(parts[p_i])
            if snp_id and p_val is not None:
                hwe_dict[snp_id] = p_val
                n_loaded += 1
            else:
                n_skipped += 1

    fh.close()
    sys.stderr.write(f"  -> {n_loaded:,} variants charges  ({n_skipped:,} lignes ignorees)\n\n")
    if n_loaded == 0:
        sys.stderr.write("ATTENTION: Aucun variant HWE charge — verifier le format du fichier.\n")
    return hwe_dict

# ─────────────────────────────────────────────────────────────
# Utilitaires
# ─────────────────────────────────────────────────────────────
def emac_from_freq(aaf: Optional[float], n: Optional[float]) -> Optional[float]:
    """EMAC = 2 * N * MAF"""
    if aaf is None or n is None:
        return None
    maf = aaf if aaf <= 0.5 else 1.0 - aaf
    return 2.0 * n * maf

# ─────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Post-filter REGENIE GWAS output (EMAC, MAF, P-value, HWE externe)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "NOTE HWE: le filtre HWE necessite --hwe-file (fichier PLINK .hwe ou .hardy).\n"
            "Sans ce fichier, utilisez --no-hwe (si HWE deja filtre en amont sur les .bed)."
        )
    )

    # I/O
    ap.add_argument("--in",  dest="inp", default="-",
                    help="Fichier REGENIE en entree (.gz ou texte)")
    ap.add_argument("--out", dest="out", default="-",
                    help="Fichier filtre en sortie (.gz ou texte)")

    # Filtres quantitatifs
    ap.add_argument("--emac-min",   type=float, default=100.0,
                    help="EMAC minimum (2*N*MAF)")
    ap.add_argument("--maf-min",    type=float, default=0.0,
                    help="MAF minimum (ex: 0.001)")
    ap.add_argument("--log10p-min", type=float, default=None,
                    help="LOG10P minimum (ex: 1.3 = p<0.05)")
    ap.add_argument("--p-max",      type=float, default=None,
                    help="P-value maximum (ex: 0.05)")

    # HWE — fichier externe (PLINK)
    hwe_grp = ap.add_mutually_exclusive_group()
    hwe_grp.add_argument(
        "--hwe-file", default=None,
        help=(
            "Fichier HWE pre-calcule par PLINK (.hwe PLINK1 ou .hardy PLINK2). "
            "Generer avec : plink --bfile mydata --hardy --out hwe_stats"
        )
    )
    hwe_grp.add_argument(
        "--no-hwe", action="store_true",
        help="Desactiver le filtre HWE (si deja filtre en amont)"
    )
    ap.add_argument("--hwe-minp", type=float, default=1e-12,
                    help="HWE p-value minimum (utilise avec --hwe-file)")

    # Colonnes REGENIE
    ap.add_argument("--id-col",     default="ID",      help="Colonne SNP ID")
    ap.add_argument("--aaf-col",    default="A1FREQ",  help="Colonne frequence allele A1")
    ap.add_argument("--n-col",      default="N",       help="Colonne taille echantillon")
    ap.add_argument("--log10p-col", default="LOG10P",  help="Colonne LOG10P")
    ap.add_argument("--sep",        default=None,
                    help="Separateur colonnes REGENIE (auto-detecte si None)")

    args = ap.parse_args()

    # ── Validation HWE ──────────────────────────────────────
    if not args.no_hwe and args.hwe_file is None:
        sys.stderr.write(
            "\nERREUR : le filtre HWE ne peut pas etre calcule depuis les summary stats "
            "REGENIE\n"
            "(les comptes de genotypes observes sont absents du fichier).\n\n"
            "Solutions :\n"
            "  A) Fournir un fichier HWE externe via --hwe-file <fichier.hwe>\n"
            "     Generer avec PLINK1 : plink  --bfile mydata --hardy --out hwe_stats\n"
            "     Generer avec PLINK2 : plink2 --pfile mydata --hardy --out hwe_stats\n\n"
            "  B) Desactiver HWE si deja filtre en amont : --no-hwe\n\n"
        )
        sys.exit(1)

    # ── Chargement HWE ──────────────────────────────────────
    hwe_dict: Dict[str, float] = {}
    if args.hwe_file:
        sys.stderr.write(f"Chargement fichier HWE : {args.hwe_file}\n")
        hwe_dict = load_hwe_file(args.hwe_file)

    # ── Ouverture fichiers ──────────────────────────────────
    inp        = open_any(args.inp)
    out_handle = open_out(args.out)

    # ── Header ──────────────────────────────────────────────
    raw_header = inp.readline().rstrip("\n")
    if not raw_header:
        sys.stderr.write("Erreur: fichier vide\n")
        return

    sep = args.sep
    if sep is None:
        sep = "\t" if "\t" in raw_header else " "

    cols    = [c.strip() for c in raw_header.split(sep)]
    col_idx = {c: i for i, c in enumerate(cols)}

    def get_idx(name):
        i = col_idx.get(name)
        if i is None:
            sys.stderr.write(f"  Colonne '{name}' non trouvee\n")
        return i

    id_i     = get_idx(args.id_col)
    aaf_i    = get_idx(args.aaf_col)
    n_i      = get_idx(args.n_col)
    log10p_i = get_idx(args.log10p_col)

    # ── Resume filtres ───────────────────────────────────────
    sys.stderr.write(f"Format  : REGENIE (sep={'TAB' if sep == chr(9) else 'ESPACE'})\n")
    sys.stderr.write(f"Colonnes: {', '.join(cols)}\n\n")
    sys.stderr.write(f"Filtres appliques:\n")
    sys.stderr.write(f"  EMAC    >= {args.emac_min}\n")
    sys.stderr.write(f"  MAF     >= {args.maf_min}\n")
    if args.no_hwe:
        sys.stderr.write(f"  HWE      : desactive (--no-hwe)\n")
    elif args.hwe_file:
        sys.stderr.write(f"  HWE p   >= {args.hwe_minp}  (source: {args.hwe_file}, "
                         f"{len(hwe_dict):,} variants)\n")
    if args.log10p_min is not None:
        sys.stderr.write(f"  LOG10P  >= {args.log10p_min}\n")
    if args.p_max is not None:
        sys.stderr.write(f"  P-value <= {args.p_max}\n")
    sys.stderr.write(f"\n")

    # ── Output header ────────────────────────────────────────
    print("\t".join(cols), file=out_handle)

    # ── Compteurs ────────────────────────────────────────────
    total           = 0
    kept            = 0
    filt_emac       = 0
    filt_maf        = 0
    filt_hwe        = 0
    filt_hwe_absent = 0
    filt_p          = 0
    filt_missing    = 0

    # ── Boucle principale ────────────────────────────────────
    for line in inp:
        line = line.rstrip("\n")
        if not line:
            continue

        total += 1
        parts = [p.strip() for p in line.split(sep)]

        snp_id = parts[id_i]  if id_i     is not None and id_i     < len(parts) else None
        aaf    = parse_float_safe(parts[aaf_i])    if aaf_i    is not None and aaf_i    < len(parts) else None
        n      = parse_float_safe(parts[n_i])      if n_i      is not None and n_i      < len(parts) else None
        log10p = parse_float_safe(parts[log10p_i]) if log10p_i is not None and log10p_i < len(parts) else None

        if aaf is None or n is None:
            filt_missing += 1
            continue

        maf  = aaf if aaf <= 0.5 else 1.0 - aaf
        emac = emac_from_freq(aaf, n)

        # ── Filtre HWE ───────────────────────────────────────
        if args.no_hwe:
            pass_hwe = True
            hwe_p    = None
        elif args.hwe_file:
            if snp_id is None or snp_id not in hwe_dict:
                filt_hwe_absent += 1
                pass_hwe = False
                hwe_p    = None
            else:
                hwe_p    = hwe_dict[snp_id]
                pass_hwe = hwe_p >= args.hwe_minp
        else:
            pass_hwe = True
            hwe_p    = None

        # ── Filtre P-value ───────────────────────────────────
        pval   = (10.0 ** (-log10p)) if log10p is not None else None
        pass_p = True
        if args.log10p_min is not None:
            pass_p = (log10p is not None) and (log10p >= args.log10p_min)
        if args.p_max is not None and pval is not None:
            pass_p = pass_p and (pval <= args.p_max)

        pass_emac = (emac is not None) and (emac >= args.emac_min)
        pass_maf  = maf >= args.maf_min

        # ── Debug 3 premiers variants ─────────────────────────
        if total <= 3:
            sys.stderr.write(f"Variant {total} ({snp_id or 'N/A'}):\n")
            sys.stderr.write(f"  AAF={aaf:.6f}  MAF={maf:.6f}  N={int(n)}  EMAC={emac:.1f}\n")
            if args.no_hwe:
                sys.stderr.write(f"  HWE=desactive\n")
            elif args.hwe_file:
                if hwe_p is not None:
                    sys.stderr.write(f"  HWE p={hwe_p:.2e} (fichier externe)  pass={pass_hwe}\n")
                else:
                    sys.stderr.write(f"  HWE=absent du fichier -> exclu\n")
            sys.stderr.write(f"  LOG10P={log10p}  P={f'{pval:.2e}' if pval else 'NA'}\n")
            sys.stderr.write(f"  -> emac={pass_emac} maf={pass_maf} hwe={pass_hwe} p={pass_p}\n\n")

        # ── Compteurs filtres ─────────────────────────────────
        if not pass_emac: filt_emac += 1
        if not pass_maf:  filt_maf  += 1
        if not pass_hwe and not (args.hwe_file and (snp_id is None or snp_id not in hwe_dict)):
            filt_hwe += 1
        if not pass_p:    filt_p    += 1

        if pass_emac and pass_maf and pass_hwe and pass_p:
            print("\t".join(parts), file=out_handle)
            kept += 1

    if out_handle is not sys.stdout:
        out_handle.close()

    # ── Stats finales ────────────────────────────────────────
    sys.stderr.write(f"\n{'='*60}\n")
    sys.stderr.write(f"Statistiques de filtrage:\n")
    sys.stderr.write(f"  Total variants              : {total:,}\n")
    sys.stderr.write(f"  Exclus (donnees manquantes) : {filt_missing:,}\n")
    sys.stderr.write(f"  Exclus (EMAC < {args.emac_min})        : {filt_emac:,}\n")
    sys.stderr.write(f"  Exclus (MAF  < {args.maf_min})         : {filt_maf:,}\n")
    if args.no_hwe:
        sys.stderr.write(f"  Exclus HWE                  : 0 (desactive)\n")
    elif args.hwe_file:
        sys.stderr.write(f"  Exclus (HWE p < {args.hwe_minp})    : {filt_hwe:,}\n")
        sys.stderr.write(f"  Exclus (absent fichier HWE) : {filt_hwe_absent:,}\n")
    if args.log10p_min is not None or args.p_max is not None:
        sys.stderr.write(f"  Exclus (seuil P)            : {filt_p:,}\n")
    sys.stderr.write(f"  Variants conserves          : {kept:,}\n")
    if total > 0:
        sys.stderr.write(f"  Taux de retention           : {100*kept/total:.2f}%\n")
    sys.stderr.write(f"{'='*60}\n")

if __name__ == "__main__":
    main()
