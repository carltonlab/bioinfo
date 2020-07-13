# bioinfo
Tools for bioinformatics made in the Carlton lab

##sequencedistance
###Finding the distance to a given feature for a given genome sequence motif
### Purpose: to be able to quickly assess the disposition of sequence motifs towards genome features such as start codons, TSS, 3' UTRs, etc

Motivation: we are interested in sequence repeats that correlate with high CO-recombination (e.g., ATTTGCCGATTTGCCG). 
These are non-randomly placed in the genome at large scale (high in arms, low in centers.)
We want to see if they are also non-randomly placed relative to exons/TSS/UTRs/etc.

There may be tools to do this already, but this project is also designed to provide a teaching opportunity for students to tie together
Emboss tools (fuzznuc), BASH, and Python.
