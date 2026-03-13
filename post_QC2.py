#!/usr/bin/env python3
"""
post_QC.py — Post-filter REGENIE GWAS output
Filtres: EMAC, HWE, MAF, P-value

Format REGENIE attendu (espace ou tab séparé):
  CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA

Usage:
  python3 post_QC.py --in gwas.regenie.gz --out gwas_filtered.tsv.gz
  python3 post_QC.py --in gwas.regenie.gz --out gwas_filtered.tsv.gz --no-hwe
  python3 post_QC.py --in gwas.regenie.gz --out gwas_filtered.tsv.gz --emac-min 50 --maf-min 0.001
"""
import sys, gzip, argparse, math
from typing import Optional

# ─────────────────────────────────────────────────────────────
# I/O
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
    except:
        return None

# ─────────────────────────────────────────────────────────────
# HWE exact test (Wigginton et al., 2005) — mid-p
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

def hwe_from_freq(aaf: float, n: int) -> Optional[float]:
    """
    Estime HWE depuis la fréquence allélique A1 et N.
    Déduit les comptes de génotypes attendus sous HWE.
    """
    if aaf is None or n is None or n == 0:
        return None
    p   = 1.0 - aaf
    q   = aaf
    exp_hom_ref = round(p * p * n)
    exp_het     = round(2 * p * q * n)
    exp_hom_alt = n - exp_hom_ref - exp_het
    if exp_hom_ref < 0 or exp_het < 0 or exp_hom_alt < 0:
        return None
    return hwe_pvalue(exp_hom_ref, exp_het, exp_hom_alt)

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
        description="Post-filter REGENIE GWAS output (EMAC, HWE, MAF, P-value)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # I/O
    ap.add_argument("--in",         dest="inp", default="-",      help="Fichier REGENIE en entrée (.gz ou texte)")
    ap.add_argument("--out",        dest="out", default="-",      help="Fichier filtré en sortie (.gz ou texte)")

    # Filtres
    ap.add_argument("--emac-min",   type=float, default=100.0,    help="EMAC minimum (2*N*MAF)")
    ap.add_argument("--maf-min",    type=float, default=0.0,      help="MAF minimum (ex: 0.001)")
    ap.add_argument("--hwe-minp",   type=float, default=1e-12,    help="HWE p-value minimum")
    ap.add_argument("--no-hwe",     action="store_true",           help="Désactiver le filtre HWE")
    ap.add_argument("--log10p-min", type=float, default=None,     help="LOG10P minimum (ex: 1.3 = p<0.05)")
    ap.add_argument("--p-max",      type=float, default=None,     help="P-value maximum (ex: 0.05)")

    # Colonnes — noms REGENIE par défaut
    ap.add_argument("--id-col",     default="ID",      help="Colonne SNP ID")
    ap.add_argument("--aaf-col",    default="A1FREQ",  help="Colonne fréquence allèle A1")
    ap.add_argument("--n-col",      default="N",       help="Colonne taille échantillon")
    ap.add_argument("--log10p-col", default="LOG10P",  help="Colonne LOG10P")
    ap.add_argument("--sep",        default=None,      help="Séparateur (auto-détecté si None)")

    args = ap.parse_args()

    inp        = open_any(args.inp)
    out_handle = open_out(args.out)

    # ── Header ──────────────────────────────────────────────
    raw_header = inp.readline().rstrip("\n")
    if not raw_header:
        sys.stderr.write("Erreur: fichier vide\n")
        return

    # Auto-détection séparateur
    sep = args.sep
    if sep is None:
        sep = "\t" if "\t" in raw_header else " "

    cols    = [c.strip() for c in raw_header.split(sep)]
    col_idx = {c: i for i, c in enumerate(cols)}

    def get_idx(name):
        i = col_idx.get(name)
        if i is None:
            sys.stderr.write(f"⚠️  Colonne '{name}' non trouvée\n")
        return i

    id_i     = get_idx(args.id_col)
    aaf_i    = get_idx(args.aaf_col)
    n_i      = get_idx(args.n_col)
    log10p_i = get_idx(args.log10p_col)

    # ── Résumé filtres ───────────────────────────────────────
    sys.stderr.write(f"\nFormat  : REGENIE (sep={'TAB' if sep == chr(9) else 'ESPACE'})\n")
    sys.stderr.write(f"Colonnes: {', '.join(cols)}\n\n")
    sys.stderr.write(f"Filtres appliqués:\n")
    sys.stderr.write(f"  EMAC    >= {args.emac_min}\n")
    sys.stderr.write(f"  MAF     >= {args.maf_min}\n")
    if args.no_hwe:
        sys.stderr.write(f"  HWE      : désactivé (--no-hwe)\n")
    else:
        sys.stderr.write(f"  HWE p   >= {args.hwe_minp}\n")
    if args.log10p_min:
        sys.stderr.write(f"  LOG10P  >= {args.log10p_min}\n")
    if args.p_max:
        sys.stderr.write(f"  P-value <= {args.p_max}\n")
    sys.stderr.write(f"\n")

    # Output header (tab-séparé)
    print("\t".join(cols), file=out_handle)

    # ── Compteurs ────────────────────────────────────────────
    total        = 0
    kept         = 0
    filt_emac    = 0
    filt_maf     = 0
    filt_hwe     = 0
    filt_p       = 0
    filt_missing = 0

    # ── Boucle ───────────────────────────────────────────────
    for line in inp:
        line = line.rstrip("\n")
        if not line:
            continue

        total += 1
        parts = [p.strip() for p in line.split(sep)]

        aaf    = parse_float_safe(parts[aaf_i])    if aaf_i    is not None and aaf_i    < len(parts) else None
        n      = parse_float_safe(parts[n_i])      if n_i      is not None and n_i      < len(parts) else None
        log10p = parse_float_safe(parts[log10p_i]) if log10p_i is not None and log10p_i < len(parts) else None

        # Données critiques manquantes → exclure
        if aaf is None or n is None:
            filt_missing += 1
            continue

        maf  = aaf if aaf <= 0.5 else 1.0 - aaf
        emac = emac_from_freq(aaf, n)

        # HWE
        if args.no_hwe:
            hwe_p    = 1.0
            pass_hwe = True
        else:
            hwe_p    = hwe_from_freq(aaf, int(n))
            pass_hwe = (hwe_p is not None) and (hwe_p >= args.hwe_minp)

        # P-value
        pval   = (10 ** (-log10p)) if log10p is not None else None
        pass_p = True
        if args.log10p_min is not None:
            pass_p = (log10p is not None) and (log10p >= args.log10p_min)
        if args.p_max is not None and pval is not None:
            pass_p = pass_p and (pval <= args.p_max)

        pass_emac = (emac is not None) and (emac >= args.emac_min)
        pass_maf  = maf >= args.maf_min

        # Debug 3 premiers variants
        if total <= 3:
            snp = parts[id_i] if id_i is not None and id_i < len(parts) else "N/A"
            sys.stderr.write(f"Variant {total} ({snp}):\n")
            sys.stderr.write(f"  AAF={aaf:.6f}  MAF={maf:.6f}  N={int(n)}  EMAC={emac:.1f}\n")
            if args.no_hwe:
                sys.stderr.write(f"  HWE=désactivé\n")
            else:
                sys.stderr.write(f"  HWE p={hwe_p:.2e}  pass={pass_hwe}\n")
            sys.stderr.write(f"  LOG10P={log10p}  P={f'{pval:.2e}' if pval else 'NA'}\n")
            sys.stderr.write(f"  → emac={pass_emac} maf={pass_maf} hwe={pass_hwe} p={pass_p}\n\n")

        if not pass_emac: filt_emac += 1
        if not pass_maf:  filt_maf  += 1
        if not pass_hwe:  filt_hwe  += 1
        if not pass_p:    filt_p    += 1

        if pass_emac and pass_maf and pass_hwe and pass_p:
            print("\t".join(parts), file=out_handle)
            kept += 1

    if out_handle is not sys.stdout:
        out_handle.close()

    # ── Stats finales ────────────────────────────────────────
    sys.stderr.write(f"\n{'='*60}\n")
    sys.stderr.write(f"Statistiques de filtrage:\n")
    sys.stderr.write(f"  Total variants             : {total:,}\n")
    sys.stderr.write(f"  Exclus (données manquantes): {filt_missing:,}\n")
    sys.stderr.write(f"  Exclus (EMAC < {args.emac_min})       : {filt_emac:,}\n")
    sys.stderr.write(f"  Exclus (MAF  < {args.maf_min})        : {filt_maf:,}\n")
    if args.no_hwe:
        sys.stderr.write(f"  Exclus HWE                 : 0 (désactivé)\n")
    else:
        sys.stderr.write(f"  Exclus (HWE p < {args.hwe_minp})   : {filt_hwe:,}\n")
    if args.log10p_min or args.p_max:
        sys.stderr.write(f"  Exclus (seuil P)           : {filt_p:,}\n")
    sys.stderr.write(f"  Variants conservés         : {kept:,}\n")
    sys.stderr.write(f"  Taux de rétention          : {100*kept/total:.2f}%\n")
    sys.stderr.write(f"{'='*60}\n")

if __name__ == "__main__":
    main()
