#!/usr/bin/env bash
PATTERN=${1? No search term given}
DNAFILE=${2:-c_elegans_ws274.fna}
GFFFILE=${3:-c_elegans_ws274_cds.gff3}

fuzznuc -comp -pat $PATTERN -seq $DNAFILE -out temporary.txt
< temporary.txt awk 'BEGIN{OFS="\t";} /pattern:/{print $1,$2,$3}' > matches.txt 
< temporary.txt awk '/HitCount:/{print $3}' > count.txt 
#rm ~/temporary.txt
    
awk 'BEGIN{OFS="\t";} /\t\+\t/{print $1,$4,$7,$8,$9}' $GFFFILE > start-posi-ex.txt
awk 'BEGIN{OFS="\t";} /\t\-\t/{print $1,$5,$7,$8,$9}' $GFFFILE > end-posi-ex.txt

python redundancy-disposal.py start-posi-ex.txt end-posi-ex.txt
#rm ~/genomicdata/start-posi-ex.txt ~/genomicdata/end-posi-ex.txt

python cds-pat.py matches.txt count.txt c_elegans_ws274_startcodons.txt $PATTERN
#rm matches.txt count.txt
