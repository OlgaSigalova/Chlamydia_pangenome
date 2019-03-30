#!/bin/bash

#
#At the AA level, MACSE represents the stop codon by its usual
#symbol ‘‘*’’ and a codon containing a frameshift is represented by
#an extra symbol, the ‘‘!’’ (see figures below for examples).
#Meanwhile, at the nucleotide level, MACSE uses the symbol ‘‘!’’
#to represent deletions of one or two nucleotides that induce
#frameshifts and it uses no special representation for the stop codon
#
#We replace these symbols with '-' to pass the files to raxml

# sed -i replaces "in place", if suffix is passed after -i back-ups are made (sometimes, '' is required)

og=$1

sed -i 's/*/-/g;s/!/-/g' /home/bsgroup/Chlamydia/fastas/${og}_macse_*
