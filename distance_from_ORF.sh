#!/usr/bin/env bash
PATTERN=${1? No search term given}
FILE=${2:-~/genomicdata/Celegans_genome.fna}
FILE2=${3:~/genomicdata/ws274.cds.txt}

fuzznuc -comp -pat $PATTERN -seq $FILE -out ~/temporary.txt
< ~/temporary.txt awk 'BEGIN{OFS="\t";} /pattern:/{print $1,$2,$3}' > matches.txt 
< ~/temporary.txt awk '/HitCount:/{print $3}' > count.txt 
rm ~/temporary.txt
    
awk 'BEGIN{OFS="\t";} /\t\+\t/{print $1,$4,$7,$8,$9}' ~/genomicdata/ws274.cds.txt > ~/genomicdata/start-posi-ex.txt
awk 'BEGIN{OFS="\t";} /\t\-\t/{print $1,$5,$7,$8,$9}' ~/genomicdata/ws274.cds.txt > ~/genomicdata/end-posi-ex.txt

python ~/scripts/2020may-genomic-analysis/20200529redundancy-disposal.py ~/genomicdata/start-posi-ex.txt ~/genomicdata/end-posi-ex.txt
rm ~/genomicdata/start-posi-ex.txt ~/genomicdata/end-posi-ex.txt

python ~/scripts/2020-June/cds-pat.py matches.txt count.txt /Users/taka/genomicdata/ws274_elegans_startcodons.txt $PATTERN
rm matches.txt count.txt
