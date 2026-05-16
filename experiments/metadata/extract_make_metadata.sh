#!/bin/bash
# metadata/extract_make_metadata.sh
# Scan make_*.R files (in ../figures_and_tables/) and build two TSVs here:
#   map_scripts.tsv - indexed by script (script -> required results -> produces)
#   map_figures.tsv - indexed by figure/table (output -> script -> required results)

set -u

cd "$(dirname "$0")" || exit 1

SRC_DIR="../figures_and_tables"
OUT=map_scripts.tsv
MAP=map_figures.tsv

printf "script\trequired_results\tproduces\n" > "$OUT"
printf "output\tscript\trequired_results\n"   > "$MAP"

for f in "$SRC_DIR"/make_*.R; do
  [[ -f "$f" ]] || continue
  base=$(basename "$f")
  awk -v file="$base" -v map="$MAP" '
    /^## Required results[[:space:]]*:/ {
      sub(/^## Required results[[:space:]]*:[[:space:]]*/, "")
      required = $0
    }
    /^## Produces[[:space:]]*:/ {
      sub(/^## Produces[[:space:]]*:[[:space:]]*/, "")
      produces = $0
      print file "\t" required "\t" produces
      n = split(produces, items, /[[:space:]]*,[[:space:]]*/)
      for (i = 1; i <= n; i++) {
        print items[i] "\t" file "\t" required >> map
      }
      required = ""; produces = ""
    }
  ' "$f" >> "$OUT"
done

{ head -n 1 "$MAP"; tail -n +2 "$MAP" | sort -V; } > "$MAP.sorted" && mv "$MAP.sorted" "$MAP"

echo "Wrote $OUT"
column -t -s $'\t' "$OUT"
echo ""
echo "Wrote $MAP"
column -t -s $'\t' "$MAP"
