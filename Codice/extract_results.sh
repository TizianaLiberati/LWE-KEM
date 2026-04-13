#!/bin/bash

OUTPUT="final_results.csv"
TMP="tmp_results.csv"

# Header
echo "n;threads;avg_keygen_us;avg_encaps_us;avg_decaps_us;total_s;mismatches" > "$TMP"

# Extract CSV lines directly
find benchmark_runs_* -type f -name "out_*.txt" | while read -r outfile; do
    line=$(grep -E '^[0-9]+;' "$outfile" | head -n1)
    [[ -n "$line" ]] && echo "$line" >> "$TMP"
done

# Sort by n
sort -t';' -k1 -n -u "$TMP" > "$OUTPUT"
rm "$TMP"

echo "✅ CSV written to $OUTPUT"

echo ""
echo "📊 Pretty Table:"
echo ""

# Pretty print with formatting
awk -F';' '
BEGIN {
    printf "%-8s %-8s %-15s %-15s %-15s %-10s %-10s\n",
           "n", "threads", "keygen(us)", "encaps(us)", "decaps(us)", "total(s)", "mismatch"
    print "---------------------------------------------------------------------------------------------"
}
NR>1 {
    # Convert scientific notation to float
    keygen = sprintf("%.2f", $3)
    encaps = sprintf("%.2f", $4)
    decaps = sprintf("%.2f", $5)
    total  = sprintf("%.6f", $6)

    printf "%-8s %-8s %-15s %-15s %-15s %-10s %-10s\n",
           $1, $2, keygen, encaps, decaps, total, $7
}
' "$OUTPUT"
