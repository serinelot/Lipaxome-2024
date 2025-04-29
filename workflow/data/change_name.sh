for f in *.fastq.gz; do 
  mv "$f" "$(echo "$f" | sed -E 's/_S[0-9]+_R([12])_001/_\1/')"
done