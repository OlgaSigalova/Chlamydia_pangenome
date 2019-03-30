opendir (DIR,'Chlamydia/fastas/universal_noframeshift_noparalogs/');
my @file= readdir DIR;
foreach $i (@file) {
	if ($i =~ /_macse_AA.fasta/)
	{
		open (F, "Chlamydia/fastas/universal_noframeshift_noparalogs/".$i);
		undef %seq;
        print "Chlamydia/fastas/universal_noframeshift_noparalogs/".$i;
		while ($str = <F>) {
			chomp $str;
            print $str;
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
open (OUT, ">Chlamydia/fastas/concatenated_core_genes_alignments_aa.fasta");
foreach $i (keys %total_seq) {
	print OUT ">".$i."\n".$total_seq{$i}."\n";
}
