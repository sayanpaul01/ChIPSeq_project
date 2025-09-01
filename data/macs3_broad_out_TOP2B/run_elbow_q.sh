#!/usr/bin/env bash
# Find per-sample elbow using qscore (−log10 FDR) within 1..3 (FDR 0.1..0.001).
# Outputs an Excel file if pandas+openpyxl are available, otherwise CSV.
# Usage:
#   chmod +x run_elbow_q.sh
#   ./run_elbow_q.sh -d . -o "$HOME/ELBOW_Q_SUMMARY.xlsx"

set -euo pipefail

SEARCH_DIR="."
OUTPUT_XLSX="ELBOW_Q_SUMMARY.xlsx"
while getopts ":d:o:" opt; do
  case "${opt}" in
    d) SEARCH_DIR="${OPTARG}";;
    o) OUTPUT_XLSX="${OPTARG}";;
    *) echo "Usage: $0 [-d DIR] [-o OUT.xlsx]"; exit 1;;
  esac
done

export SEARCH_DIR OUTPUT_XLSX

python3 - <<'PY'
import os, sys, glob, csv, math

SEARCH_DIR  = os.environ.get("SEARCH_DIR", ".")
OUTPUT_XLSX = os.environ.get("OUTPUT_XLSX", "ELBOW_Q_SUMMARY.xlsx")

# Try Excel deps; fall back to CSV
have_pandas = have_openpyxl = False
try:
    import pandas as pd
    have_pandas = True
    try:
        import openpyxl  # noqa
        have_openpyxl = True
    except Exception:
        have_openpyxl = False
except Exception:
    have_pandas = False

def read_q_np_lp(fn):
    """Return lists (qscore, npeaks, lpeaks). Skip file if needed cols are missing."""
    qs, nps, lps = [], [], []
    with open(fn, newline='') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        if not rdr.fieldnames: return qs, nps, lps
        cols = [c.strip() for c in rdr.fieldnames]
        if "qscore" not in cols or "npeaks" not in cols or "lpeaks" not in cols:
            return qs, nps, lps
        for r in rdr:
            try:
                q  = float(r["qscore"])
                np = float(r["npeaks"])
                lp = float(r["lpeaks"])
            except Exception:
                continue
            qs.append(q); nps.append(np); lps.append(lp)
    return qs, nps, lps

def elbow_by_secant(xs, ys, lo, hi):
    """
    xs (qscore), ys (lpeaks).
    1) keep lo <= x <= hi
    2) sort by x
    3) scale x,y to [0,1] within the window
    4) find point with max perpendicular distance to secant between endpoints
    returns x_at_elbow or None
    """
    pts = [(x, y) for x, y in zip(xs, ys) if lo <= x <= hi]
    if len(pts) < 3:
        return None
    pts.sort(key=lambda t: t[0])
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    dx = (xmax - xmin) or 1e-12
    dy = (ymax - ymin) or 1e-12

    # Normalize
    xn = [(x - xmin)/dx for x in xs]
    yn = [(y - ymin)/dy for y in ys]

    # Secant endpoints in normalized space
    x0, y0 = xn[0], yn[0]
    x1, y1 = xn[-1], yn[-1]
    denom = math.hypot(y1 - y0, x1 - x0) or 1e-12

    # Perpendicular distance of each point to the secant
    def dist(i):
        x, y = xn[i], yn[i]
        # |(y1 - y0)x - (x1 - x0)y + x1*y0 - y1*x0| / sqrt((y1-y0)^2 + (x1-x0)^2)
        return abs((y1 - y0)*x - (x1 - x0)*y + x1*y0 - y1*x0)/denom

    best_i, best_d = 0, -1.0
    for i in range(len(xs)):
        d = dist(i)
        if d > best_d:
            best_i, best_d = i, d
    return xs[best_i]  # return the original (un-normalized) qscore at the elbow

def nearest_row(xs, nps, lps, target):
    """Return (x*, np*, lp*) for the row closest to target."""
    j = min(range(len(xs)), key=lambda i: abs(xs[i] - target))
    return xs[j], nps[j], lps[j]

MIN_Q, MAX_Q = 1.0, 3.0  # qscore window == −log10(FDR) in [1,3]

files = sorted(glob.glob(os.path.join(SEARCH_DIR, "**", "*_cutoff_analysis.txt"), recursive=True))
if not files:
    sys.exit(f"No *_cutoff_analysis.txt files found under: {SEARCH_DIR}")

rows = []
for fn in files:
    sample = os.path.basename(fn).replace("_cutoff_analysis.txt", "")
    qs, nps, lps = read_q_np_lp(fn)
    if not qs:
        continue

    elbow = elbow_by_secant(qs, lps, MIN_Q, MAX_Q)
    if elbow is None:
        # fallback: if nothing in [1,3], try all q >=1; else max q
        eligible = [(q, lp) for q, lp in zip(qs, lps) if q >= 1.0]
        if len(eligible) >= 3:
            ss, vv = zip(*eligible)
            elbow = elbow_by_secant(list(ss), list(vv), min(ss), max(ss))
        if elbow is None:
            elbow = max(qs)

    q_star, np_star, lp_star = nearest_row(qs, nps, lps, elbow)
    np_out = int(round(np_star)) if abs(np_star - round(np_star)) < 1e-9 else float(np_star)

    rows.append({
        "Sample": sample,
        "Elbow_qscore": round(float(q_star), 6),
        "NPeaks_at_elbow": np_out,
        "LPeaks_at_elbow": float(lp_star),
        # Optional: raw FDR if you want to keep an eye on it
        "FDR_at_elbow": 10**(-float(q_star)),
    })

if not rows:
    sys.exit("No valid rows computed.")

# Write output
if have_pandas and have_openpyxl:
    import pandas as pd  # type: ignore
    df = pd.DataFrame(rows).sort_values("Sample")
    try:
        df.to_excel(OUTPUT_XLSX, index=False)
        print(f"Wrote {len(df)} rows to {OUTPUT_XLSX}")
    except PermissionError:
        home_out = os.path.join(os.path.expanduser("~"), os.path.basename(OUTPUT_XLSX))
        df.to_excel(home_out, index=False)
        print(f"Permission denied for '{OUTPUT_XLSX}'. Wrote to {home_out}")
else:
    out_csv = os.path.splitext(OUTPUT_XLSX)[0] + ".csv"
    try:
        with open(out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Sample","Elbow_qscore","NPeaks_at_elbow","LPeaks_at_elbow","FDR_at_elbow"])
            for r in sorted(rows, key=lambda r: r["Sample"]):
                w.writerow([r["Sample"], r["Elbow_qscore"], r["NPeaks_at_elbow"], r["LPeaks_at_elbow"], r["FDR_at_elbow"]])
        print(f"pandas/openpyxl not found — wrote CSV instead: {out_csv}")
        print("Install Excel deps with:  python3 -m pip install --user pandas openpyxl")
    except PermissionError:
        home_out = os.path.join(os.path.expanduser("~"), os.path.basename(out_csv))
        with open(home_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Sample","Elbow_qscore","NPeaks_at_elbow","LPeaks_at_elbow","FDR_at_elbow"])
            for r in sorted(rows, key=lambda r: r["Sample"]):
                w.writerow([r["Sample"], r["Elbow_qscore"], r["NPeaks_at_elbow"], r["LPeaks_at_elbow"], r["FDR_at_elbow"]])
        print(f"Permission denied for '{out_csv}'. Wrote CSV to {home_out}")
PY
