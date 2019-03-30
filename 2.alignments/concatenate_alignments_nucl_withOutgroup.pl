opendir (DIR,'Chlamydia/fastas/universal_noframeshift_noparalogs_withOutgroup/');
my @file= readdir DIR;
foreach $i (@file) {
	if ($i =~ /_macse_NT.fasta/)
	{
		open (F, 'Chlamydia/fastas/universal_noframeshift_noparalogs_withOutgroup/'.$i);
		undef %seq;
		while ($str = <F>) {
			chomp $str;
			if ($str =~ />(.*)$/){
				$curr_name = substr($1,0,4);
			} else {
				$seq{$curr_name} .= $str; 
			}
		}
		foreach $i (keys %seq) {
			$total_seq{$i} .= $seq{$i};
		}
		close F;
	}
}
closedir DIR;
open (OUT, ">Chlamydia/fastas/concatenated_core_genes_alignments_nucl_withOutgroup.fasta");
foreach $i (keys %total_seq) {
	print OUT ">".$i."\n".$total_seq{$i}."\n";
}
