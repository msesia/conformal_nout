#!/bin/bash
# metadata/extract_submit_metadata.sh
# Scan submit_*.sh files (in the parent directory) and build two TSVs here:
#   map_submit.tsv   - indexed by script (script -> produces -> output_dir)
#   map_results.tsv  - indexed by result group (result -> script -> output_dir)

set -u

cd "$(dirname "$0")" || exit 1

SRC_DIR=".."
OUT=map_submit.tsv
MAP=map_results.tsv

printf "script\tproduces\toutput_dir\n"  > "$OUT"
printf "result\tscript\toutput_dir\n"    > "$MAP"

for f in "$SRC_DIR"/submit_*.sh; do
  [[ -f "$f" ]] || continue
  base=$(basename "$f")
  awk -v file="$base" -v map="$MAP" '
    /^## Produces results[[:space:]]*:/ {
      sub(/^## Produces results[[:space:]]*:[[:space:]]*/, "")
      produces = $0
      in_produces = 1
      next
    }
    /^## Output directory[[:space:]]*:/ {
      sub(/^## Output directory[[:space:]]*:[[:space:]]*/, "")
      output_dir = $0
      in_produces = 0

      print file "\t" produces "\t" output_dir
      n = split(produces, items, /[[:space:]]*,[[:space:]]*/)
      for (i = 1; i <= n; i++) {
        if (items[i] == "") continue
        dir = output_dir
        gsub(/<FIGURE>/, items[i], dir)
        gsub(/<TABLE>/,  items[i], dir)
        print items[i] "\t" file "\t" dir >> map
      }
      produces = ""; output_dir = ""
      next
    }
    in_produces && /^##/ && !/^## -+/ {
      line = $0
      sub(/^##[[:space:]]*/, "", line)
      produces = produces " " line
      next
    }
    /^## -+/ { in_produces = 0 }
  ' "$f" >> "$OUT"
done

{ head -n 1 "$MAP"; tail -n +2 "$MAP" | sort -V; } > "$MAP.sorted" && mv "$MAP.sorted" "$MAP"

echo "Wrote $OUT"
column -t -s $'\t' "$OUT"
echo ""
echo "Wrote $MAP"
column -t -s $'\t' "$MAP"
