
cd Chlamydia/raxml/
cp /home/bsgroup/Chlamydia/fastas/concatenated_core_genes_alignments_nucl.fasta ./

/home/bsgroup/soft/standard-RAxML-8.2.3/raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -T 16 -s concatenated_core_genes_alignments_nucl.fasta -n core_genes_phylogeny

rm concatenated_core_genes_alignments_nucl.fasta
