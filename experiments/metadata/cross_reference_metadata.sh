#!/bin/bash
# metadata/cross_reference.sh
# Join the four metadata TSVs in this directory into a single end-to-end map:
#
#   figure/table  <-  required result  <-  make script
#                                      ->  submit script  ->  result output dir
#
# Inputs (in this directory):
#   map_submit.tsv     script -> produces -> output_dir
#   map_results.tsv    result -> script -> output_dir
#   map_scripts.tsv    make script -> required results -> produces
#   map_figures.tsv    figure/table -> make script -> required result
#
# Outputs (in this directory):
#   map_full.tsv  Paper item | Submit script | Results directory | Make script | Result key
#   map_full.md   same content as a Markdown table (for README / ACC form)

set -u

# Always run from this script's directory so paths resolve consistently.
cd "$(dirname "$0")" || exit 1

MAP_SUBMIT="map_submit.tsv"
MAP_RESULTS="map_results.tsv"
MAP_SCRIPTS="map_scripts.tsv"
MAP_FIGURES="map_figures.tsv"

OUT_TSV="map_full.tsv"
OUT_MD="map_full.md"

for f in "$MAP_SUBMIT" "$MAP_RESULTS" "$MAP_SCRIPTS" "$MAP_FIGURES"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: $f not found in $(pwd)" >&2
    exit 1
  fi
done

awk -F'\t' -v OFS='\t' '
  NR == FNR {
    if (FNR == 1) next
    submit_script[$1] = $2
    result_dir[$1]    = $3
    next
  }
  FNR == 1 {
    print "Paper item", "Submit script", "Results directory", "Make script", "Result key"
    next
  }
  {
    output          = $1
    make_script     = $2
    required_result = $3
    s = (required_result in submit_script) ? submit_script[required_result] : "(none)"
    d = (required_result in result_dir)    ? result_dir[required_result]    : "(none)"
    print output, s, d, make_script, required_result
  }
' "$MAP_RESULTS" "$MAP_FIGURES" > "$OUT_TSV"

{ head -n 1 "$OUT_TSV"; tail -n +2 "$OUT_TSV" | sort -V; } > "$OUT_TSV.sorted" \
  && mv "$OUT_TSV.sorted" "$OUT_TSV"

awk -F'\t' '
  NR == 1 {
    n = NF
    line = "|"; sep = "|"
    for (i = 1; i <= n; i++) {
      line = line " " $i " |"
      sep  = sep  "---|"
    }
    print line
    print sep
    next
  }
  {
    line = "|"
    for (i = 1; i <= NF; i++) line = line " " $i " |"
    print line
  }
' "$OUT_TSV" > "$OUT_MD"

OUT_TEX="map_full.tex"

awk -F'\t' '
  function escape(s) {
    gsub(/_/, "\\_", s)
    return s
  }
  NR == 1 {
    n = NF
    spec = ""
    for (i = 1; i <= n; i++) spec = spec "l"
    print "\\begin{tabular}{" spec "}"
    print "\\toprule"
    line = ""
    for (i = 1; i <= n; i++) {
      line = line escape($i)
      if (i < n) line = line " & "
    }
    print line " \\\\"
    print "\\midrule"
    next
  }
  {
    line = ""
    for (i = 1; i <= NF; i++) {
      cell = escape($i)
      if (i == 1) {
        line = line cell                   # Paper item: plain text
      } else {
        line = line "\\texttt{" cell "}"   # other columns: typewriter
      }
      if (i < NF) line = line " & "
    }
    print line " \\\\"
  }
  END {
    print "\\bottomrule"
    print "\\end{tabular}"
  }
' "$OUT_TSV" > "$OUT_TEX"

echo "Wrote $OUT_TEX (LaTeX table for response letter / appendix)"


echo "Wrote $OUT_TSV"
echo ""
column -t -s $'\t' "$OUT_TSV"
echo ""
echo "Wrote $OUT_MD (Markdown table for README / ACC form)"

# --- Diagnostics ---
echo ""
echo "============================================================"
echo "Diagnostics"
echo "============================================================"

awk -F'\t' 'NR==FNR { if (FNR>1) used[$3]=1; next } FNR>1 && !($1 in used) { print "  - " $1 " (produced by " $2 ", not consumed)" }' \
  "$MAP_FIGURES" "$MAP_RESULTS" \
  | { mapfile -t orphans; if [[ ${#orphans[@]} -gt 0 ]]; then
        echo "Unused results (produced but no make_*.R requires them):"
        printf '%s\n' "${orphans[@]}"
      else
        echo "All produced results are used by at least one make_*.R."
      fi; }

awk -F'\t' 'NR==FNR { if (FNR>1) prod[$1]=1; next } FNR>1 && !($3 in prod) { print "  - " $3 " (required by " $2 " for " $1 ")" }' \
  "$MAP_RESULTS" "$MAP_FIGURES" \
  | { mapfile -t missing; if [[ ${#missing[@]} -gt 0 ]]; then
        echo ""
        echo "Missing results (required but not produced by any submit_*.sh):"
        printf '%s\n' "${missing[@]}"
      fi; }
