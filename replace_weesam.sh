 1811  zero=$(bedtools genomecov -bga -ibam mapped_to_HCV_database_sorted.bam | awk '$4==0 {bpCountZero+=($3-$2)} {print bpCountZero}' | tail -1)
 1812  nonzero=$(bedtools genomecov -bga -ibam mapped_to_HCV_database_sorted.bam | awk '$4>0 {bpCountNonZero+=($3-$2)} {print bpCountNonZero}' | tail -1)
 1813  percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero))*100")
 1814  echo $percent

