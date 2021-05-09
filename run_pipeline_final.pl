#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;
use Bio::Perl;
#use Config::File;
use File::Basename;

my ($out_dir, $genome_id, $cov, $exp, $start, $test, $help, $all_rsem, $blast, $res_dir,$exp_cutoff, $int_id, $mat_id);
$genome_id = "GEN";
my $cov_cutoff = 0.2;
$exp_cutoff = 1;
my $exp_type = "best";
my $exp_prog = "rsem_tpm";
$res_dir = cwd();
GetOptions("m|matrix:s"=>\$mat_id, "r|result:s"=>\$res_dir, "c|cov:s"=>\$cov, "i|id:s"=>\$int_id, "e|exp:s"=>\$exp_cutoff, "g|genome:s"=>\$genome_id,  "o|out_dir:s"=>\$out_dir, "b|blast_dir:s"=>\$blast, "h|?|help"=>\$help);

if ($help || !$blast)
{
    print STDERR "Runs Post-Redundancy Removal Trimming on a Trinity Sample\n\n";
    print STDERR "-----------------USAGE--------------\n";
    print STDERR " -b database \n";
    print STDERR " -r results directory\n";
    print STDERR " -o output directory. Default is the current directory\n";
    print STDERR " -g genome id. Default = GEN\n";
    print STDERR " -i intermediate id. Default = GEN\n";
    print STDERR " -c gene coverage cutoff for good. Default = 0.2\n";
    print STDERR " -m expression matrix ID. Default = intermediate ID.\n";
    print STDERR " -e expression cutoff for good. Default = 1 FPKM\n";
    
    exit();
}


my $config_file = dirname($ARGV[0]) . "/trtap.ini";
my $config_hash ;
open(my $conf_file, "<", $config_file);
while (<$conf_file>) { chomp; # no newline 
        s/#.*//; # no comments 
        s/^\s+//; # no leading 
        s/\s+$//; # no trailing white next 
        next unless length; # anything left? 
        my ($var, $value) = split(/\s*=\s*/, $_, 2);
        $config_hash->{$var} = $value;
}


#read in gene_spreadsheet
my %list;
my %allgo;
my %overlap;
if (!$mat_id)
{
	$mat_id = $int_id;
}
#id, score, allele, match_genome,  match_id, match_len, len_cov, e_value, start_pos, start, dir, posschi, second_best, bit_score, translate_len
my @head = ("id", "score", "allele", "match_genome",  "match_id", "match_len", "len_cov", "e_value", "start_pos", "orig_len", "start", "dir", "posschi", "second_best_allele", "bit_score", "translate_len");
if (-e "$res_dir/$genome_id" .".gene_info.txt")
{
	open(my $fi, "<", "$res_dir/$genome_id" .".gene_info.txt"); <$fi>;
	while (my $v = <$fi>)
	{
		$v =~ /([^\n\r]+)/; my @n = split "\t", $1;
		for (my $i = 0; $i < scalar(@n); $i++)
		{
			$list{$n[0]}->{$head[$i]} = $n[$i];
		}
	}
	close($fi);
}
else
{
	warn("Cannot identify gene info file... quitting");
}
foreach my $a (keys(%list))
{
	if ($list{$a}->{match_genome} && $list{$a}->{match_id})
	{
		
		$overlap{$list{$a}->{match_genome} ."-" . $list{$a}->{match_id}}->{list} .= $a .";";
		if ($list{$a}->{score} eq "BEST")
		{
			$overlap{$list{$a}->{match_genome} ."-" . $list{$a}->{match_id}}->{best} = $a;
		}
	}
}
#check for rsem_matrix/rsem and bowtie
if (-e "$res_dir/$mat_id" .".rsem_matrix.tpm.txt")
{
	%list = %{read_all_rsem(\%list, "$res_dir/$mat_id" .".rsem_matrix.tpm.txt", "rsem_tpm")};
	if (-e "$res_dir/$mat_id" .".bowtie_matrix.unamb_goodcnts.txt")
	{
		%list = %{read_all_rsem(\%list, "$res_dir/$mat_id" .".bowtie_matrix.unamb_goodcnts.txt", "botwtie_cnts")};   
	}
}
else
{
	warn("Cannot find expression files $res_dir/$mat_id" .".rsem_matrix.tpm.txt");
}
print STDERR "Getting Blast File Data\n";
my $ls_sys = "ls $res_dir/$int_id" . "*blast*";
my $ls = `$ls_sys`;
while ($ls =~ /([^\n\r]+)/g)
{
	my $f = $1;
	if (-e $f)
	{
		print STDERR "$f\n";
		open(my $fb, "<", $f) or die();
		while (my $v = <$fb>)
		{
			$v =~ /([^\n\r]+)/; my @n = split "\t", $1;
			$n[0] =~ /(\S+)_i(\d+)/; my $id = $1; my $allele = $2;
			if ($list{$id}->{match_id} eq $n[1] && $list{$id}->{allele} == $allele)
			{
				#die ($n[0] . " " . $n[8] . " " . $n[9]);
				if ($n[8] < $n[9])
				{
					$list{$id}->{bl_st} = $n[8]; $list{$id}->{bl_end} = $n[9]; 
				} else {  $list{$id}->{bl_st} = $n[9]; $list{$id}->{bl_end} = $n[8]; }
				$list{$id}->{bl_range} = min($n[8], $n[9]) . "-" . max($n[8], $n[9]);
			}
		}
	}
}
print STDERR "Getting Overlap Data\n";
foreach my $a (keys(%overlap))
{
	my @lst = split ";", $overlap{$a}->{list}; my $bst = $overlap{$a}->{best}; my $b_s = $list{$bst}->{bl_st}; my $b_e = $list{$bst}->{bl_end};
	if ($bst)
	{	
		foreach my $b (@lst)
		{
			#die($b. " " . $list{$b}->{bl_end} . " " .  $list{$b}->{bl_st});
			if ($list{$b}->{bl_end}  && $list{$b}->{bl_st})
			{		
			if ($list{$b}->{bl_st} > $b_e || $list{$b}->{bl_end} < $b_s)
			{
				$list{$b}->{overlap} = "NO OVERLAP";
			} 
			else
			{
				$list{$b}->{overlap} = 100 * (min($list{$b}->{bl_end}, $b_e) - max($list{$b}->{bl_st}, $b_s)) / ($list{$b}->{bl_end} - $list{$b}->{bl_st});
			}
		
			}
		
		}
	}
}
if (-e "$res_dir/$genome_id" .".hmm.txt")
{
	my %all_pfam;
	open(my $f_all, "<", "$res_dir/$genome_id". ".hmm.txt");
	while(my $v = <$f_all>)
	{
		$v =~ /([^\n\r]+)/; my @n = split " ", $1;
		if ($n[1] =~ /PF(\d+)/)
		{
			$all_pfam{"PF$1"} = 1;
		}
	}
	my %go_pfam;
	open(my $f_go, "<", "$blast/PFAM_gene_ontology.txt");
	while(my $v = <$f_go>)
	{
		$v =~ /([^\n\r]+)/; my @n = split "\t", $1;
		if ($all_pfam{$n[0]})
		{
			$go_pfam{$n[0]}->{$n[1]} = 1;
		}
	}
	my %def_pfam;
	open(my $f_go, "<", "$blast/pfamA.txt");
	while(my $v = <$f_go>)
	{
		$v =~ /([^\n\r]+)/; my @n = split "\t", $1;
		if ($all_pfam{$n[0]})
		{
			$def_pfam{$n[0]}->{$n[3]} = 1;
		}
	}
	open(my $f_all, "<", "$res_dir/$genome_id". ".hmm.txt") or warn("Cannot open HMM file  " ."$res_dir/$int_id". ".hmm.txt");
	while(my $v = <$f_all>)
	{
		$v =~ /([^\n\r]+)/; my @n = split " ", $1;
		if ($n[1] =~ /PF(\d+)/)
		{	
			my $num = $1;
			foreach my $a (keys(%{$go_pfam{"PF$num"}}))
			{
				my $tmp = "$a;";
				if ($list{$n[2]}->{pfam_go} !~ /$tmp/)
				{
					$allgo{$n[2]}->{$tmp}++;
					$list{$n[2]}->{pfam_go} .= "$tmp";
				}
			}
			foreach my $a (keys(%{$def_pfam{"PF$num"}}))
			{
				my $tmp = "$a;";
				
				$list{$n[2]}->{pfam_def} .= "$tmp";
				
			}
			my $tmp = $n[1] . ";";
			if ($list{$n[2]}->{pfam_id} !~ /$tmp/)
			{
				$list{$n[2]}->{pfam_id} .= "$tmp";
			}
			
		}
	}
}
else
{
        warn("Cannot find PFAM file..");
}
if (-e "$res_dir/".$genome_id. "_v_SwissProt.blastx.txt")
{
	my %sw_tmp;
	my %sw_tax;
	open(my $f_sw, "<", "$blast/uniprot_sprot.ids") or die();
	while(my $v = <$f_sw>)
	{
		$v =~ /([^\r\n]+)/; my @n = split "\t", $1;
		while($n[1] =~ /([^\s;]+)/g)
		{
			$sw_tmp{$1}->{tax} = $n[3];
			$sw_tax{$n[3]}->{id} = 1;
			$sw_tmp{$1}->{name} = $n[2];
			$sw_tmp{$1}->{GO} = $n[4];
		}
	}
	close($f_sw);
	my $curr;
	open(my $sw_dat, "<",  "$blast/uniprot_sprot.dat");
	while(my $v = <$sw_dat>)
	{
		$v =~ /\A(\w{2})(\s+)([^\n\r]+)/; my $id1 = $1; my $def = $3;
		if ($id1 eq "\/\/") { $curr= ""; }
		if ($id1 eq "AC") { $def =~ /\A([^\s\;]+)/; $curr = $1; }
		if ($id1 eq "OC" && $curr && $sw_tax{$curr}->{id}) { $sw_tax{$curr}->{full} = $def; }
	}
	close($sw_dat);
	open(my $sw_tx, "<",  "$blast/taxonomy.xml");
	while (my $v = <$sw_tx>)
	{
		if ($v =~ /\A\</)
		{
			$v =~ /taxId\=\"(\d+)\"/; my $tid = $1;
			if ($sw_tax{$tid}->{id})
			{
				if ($v =~ /scientificName\=\"([^\"]+)\"/) {  $sw_tax{$curr}->{tax_name} = $1; }
				if ($v =~ /taxonomicDivision\=\"([^\"]+)\"/) { $sw_tax{$curr}->{tax_div} = $1; }
			}
		}
	}
	open(my $bl_sw, "<",  "$res_dir/".$genome_id. "_v_SwissProt.blastx.txt");
	while(my $v = <$bl_sw>)
	{
	 	$v =~ /([^\r\n]+)/; my @n = split "\t", $1;
                if ($sw_tmp{$n[1]})
                {
			$list{$n[0]}->{sw_id} = $n[1];
			$list{$n[0]}->{sw_go} = $sw_tmp{$n[1]}->{GO};
			while ($sw_tmp{$n[1]}->{GO} =~ /([A-Za-z0-9:]+)/g) { $allgo{$n[0]}->{$1.";"}++; }
			$list{$n[0]}->{sw_def} = $sw_tmp{$n[1]}->{name};
			$list{$n[0]}->{sw_tax} = $sw_tmp{$n[1]}->{tax};
			$list{$n[0]}->{sw_tax_full} = $sw_tax{$n[1]}->{full};
			$list{$n[0]}->{sw_tax_name} = $sw_tax{$n[1]}->{tax_name};
			$list{$n[0]}->{sw_tax_div} = $sw_tax{$n[1]}->{tax_div};
		}
	}
} 
my $busco_id = $genome_id ."_BUSCO";
if (!-e "$res_dir/run_$busco_id" . "/full_table_$busco_id" . ".tsv") { $busco_id = $genome_id; }
if (!-e "$res_dir/run_$busco_id" . "/full_table_$busco_id" . ".tsv") { $busco_id = $mat_id; }
if (!-e "$res_dir/run_$busco_id" . "/full_table_$busco_id" . ".tsv") { $busco_id = $mat_id . "_BUSCO"; }

if (-e "$res_dir/run_$busco_id" . "/full_table_$busco_id" . ".tsv")
{
	open(my $fb, "<", "$res_dir/run_$busco_id" . "/full_table_$busco_id" . ".tsv");
	while(my $v = <$fb>)
	{
		$v =~ /([^\n\r]+)/; my @n = split "\t", $1;
		if ($list{$n[2]})
		{
			$list{$n[2]}->{BUSCO} = $n[0] . ":" . $n[1];
		}
	}	
}
else
{
	warn("Cannot find BUSCO... $busco_id\n");
}
if (-e "$res_dir/$int_id" ."_v_Fly.blastx.out")
{
	my %fly_map;
	open(my $fm, "<", "$blast/fbgn_fbtr_fbpp_fb_2019_06.tsv");
	while(my $v = <$fm>)
	{
		$v =~ /([^\n\r]+)/g; my @n = split "\t", $1;
		if ($n[2] =~ /FBpp/)
		{
			$fly_map{$n[2]} = $n[0];
		}
	}
	my %go_map;
	open(my $fa, "<", "$blast/gene_association.fb");
	while(my $v = <$fa>)
	{
			$v =~ /([^\n\r]+)/g; my @n = split "\t", $1;
			$go_map{$n[1]}->{$n[4]} = 1;
	}
	my %sum_map;
	open(my $fa, "<", "$blast/automated_gene_summaries.tsv");
	while(my $v = <$fa>)
	{
			$v =~ /([^\n\r]+)/g; my @n = split "\t", $1;
			$sum_map{$n[0]} = $n[1];
	}
	
	my %map;
	open(my $fb, "<", "$res_dir/$int_id" ."_v_Fly.blastx.out");
	while(my $v = <$fb>)
	{
		$v =~ /([^\n\r]+)/g; my @n = split "\t", $1;
		$n[0] =~ /(\S+)_i(\d+)/;
		my $allele = $2;
		my $gene = $1;
		if ($list{$gene}->{allele} eq $allele)
		{
			my $b = $fly_map{$n[1]};
			$list{$gene}->{fly_summary} = $sum_map{$b};
		
			foreach my $a (keys(%{$go_map{$b}}))
			{
				my $tmp = "$a;";
				$allgo{$gene}->{$tmp}++;
				if ($list{$gene}->{fly_go} !~ /$tmp/)
				{
					$list{$gene}->{fly_go} .= "$tmp";
				}
			}
		}
	}
}
else
{
	warn("Cannot find FLy file..");
}
my %blcl;
warn "Getting BlastClust Data...\n";
if (-e "$res_dir/$genome_id".".blcl.95.95.out")
{
	open(my $fbc, "<", "$res_dir/$genome_id".".blcl.95.95.out");
	my $c = 1; 
	while (my $v = <$fbc>)
	{
		$v =~ /([^\n\r]+)/; my @n = split " ", $1;
		foreach my $a (@n)
		{
			$list{$a}->{blcl_num} = $c; $list{$a}->{blcl_cnt} = scalar(@n); $list{$a}->{blcl_best} = $n[0];	
		}
		$c++;
	}
}
else
{
	if (-e "$res_dir/$genome_id".".blcl.90.95.out")
	{
        open(my $fbc, "<", "$res_dir/$genome_id".".blcl.90.95.out");
        my $c = 1;
        while (my $v = <$fbc>)
        {
                $v =~ /([^\n\r]+)/; my @n = split " ", $1;
                foreach my $a (@n)
                {
                        $list{$a}->{blcl_num} = $c; $list{$a}->{blcl_cnt} = scalar(@n); $list{$a}->{blcl_best} = $n[0];
                }
                $c++;
        }
	}
}
print STDERR "Getting all GO values...\n";
my %comb_go;
if (%allgo)
{
	open(my $fgo, ">", "$res_dir/$int_id". "_all_go.txt");
	foreach my $a (keys(%allgo))
	{
		foreach my $b (keys(%{$allgo{$a}}))
		{
			$list{$a}->{all_go} .= "$b";
			$b =~ s/;//;
			print {$fgo} "$a\t$b\n";
			#print {$fgo} "ALLGO\t$a\t$a\t0\t$b\t", $allgo{$a}->{$b.";"}, "\n";
		}
	
	}
	close($fgo);
	my $sys = "cd /bigdata/hayashilab/shared/ayoublab/toby/go-perl\nperl ./goslimviewer_standalone.pl -i $res_dir/$int_id". "_all_go.txt -s generic -o $res_dir/$int_id". "_all_go_slim";
	system($sys);
	if (-e "$res_dir/$int_id". "_all_go_slim.generic.s2p2g.txt")
	{
		open(my $fgo, "<", "$res_dir/$int_id". "_all_go_slim.generic.s2p2g.txt");
		while (my $v = <$fgo>)
		{
			$v =~ /([^\n\r]+)/; my @n = split "\t", $1;
			my $b = $n[1] . ";"; if ($list{$n[3]}->{go_slim} !~ /$b/) { $list{$n[3]}->{go_slim} .= "$b"; }
		}
	}
}	
foreach my $a (keys(%list))
{
	if ($list{$a}->{score} eq "NOT_BEST")
	{
		if ($list{$a}->{$exp_type ."_" . $exp_prog} < $exp_cutoff)
		{
			$list{$a}->{score} = "LOW_EXP";
		}
		else
		{
			if ($list{$a}->{"len_cov"} < $cov_cutoff)
			{
				$list{$a}->{score} = "LOW_COV";
			}
			else
			{
				$list{$a}->{score} = "GOOD";
			}
		}
	}
	if ($list{$a}->{score} eq "LONGORF")
        {
                if ($list{$a}->{$exp_type ."_" . $exp_prog} < $exp_cutoff)
                {
                        $list{$a}->{score} = "LOW_EXP_LONGORF";
                }
	}
}
my $out_id = $genome_id;
if ($out_id =~ /INT/) { $out_id =~ s/INT/FIN/; } else { $out_it .= "_FIN"; }

open(my $fo, ">", $out . "/" . $out_id . ".gene_info.out");
my %tot;
my %type;
my %uni;
my @h = ("id", "score", "allele", "match_genome",  "match_id", "match_len", "len_cov", "e_value", "start_pos", "start", "dir", "posschi", "second_best_allele", "bit_score", "orig_len", "translate_len", "fly_go", "pfam_def", "pfam_go", "sw_def", "sw_tax", "sw_tax_full", "sw_tax_name", "sw_tax_div", "sw_go","BUSCO", "best_rsem_tpm", "best_rsem_tpm_id", "tot_rsem_tpm", "best_median_rsem_tpm", "best_median_lib_rsem_tpm", "best_botwtie_cnts", "best_botwtie_cnts_id", "tot_botwtie_cnts", "all_go", "go_slim", "blcl_num", "blcl_cnt", "blcl_best", "bl_range", "overlap");
my @hread = ("Trinity ID","Score","Trinity Allele Selected","DB Genome with Best Match","ID of the Best Matching DB Gene","Full Protein Length of the Best Match GENE","Length of the Coverage","E-Value of BlastX Martch","Frame Start Position (1,2 or 3)","Direction of the Hit","Is a Possible Chimeric","Second Best Allele","Bit Score of BlastX Hit","Length of Trinity Sequence","Length of Translated Trinity","GO Terms from Drosophila Match","Name of PFAM ID","GO Terms from PFAM","SwissProt Definition","NCBI Taxonomy ID from SwissProt Match","Taxonomy String from SwissProt", "Taxonony Name from SwissProt", "Taxonomy Division from SwissProt", "GO Terms from Swiss Prot Match","BUSCO ID Match","Highest RSEM TPM","Library with the Highest RSEM TPM","Total RSEM TPM over all Libraries","Highest Tissue Median RSEM TPM","Tissue with the Highest Median RSEM TPM","Highest Bowtie2 Count","Library with the Highest Bowtie2 Count","Total Bowtie2 Count over all Libraries","Combined GO Terms from Drosophila, SwissProt and PFAM","Generic GO Slim terms from Combined GO Terms", "BlastClust Cluster ID", "Number of Contigs in BlastClust", "Best Contig in BlastClust Cluster", "Blast Range", "% overlap with BEST"); 
foreach my $h (@h) { print {$fo} "$h\t"; } print {$fo} "\n";
foreach my $a (keys(%list))
{
	if ($list{$a}->{id} && $list{$a}->{id} ne "NA")
	{
		my $out1;
		foreach my $h (@h) { 
			if ($list{$a}->{$h}) { $out1 .=  $list{$a}->{$h} . "\t"; 
				if (!$uni{$list{$a}->{score}}) { $uni{$list{$a}->{score}} = 1; }
				if ($list{$a}->{$h} !~ /([^0-9\-\.]+)/) { $type{$h}->{$list{$a}->{score}} += $list{$a}->{$h}; $tot{$h} += $list{$a}->{$h}; 
				} else { 
					if ($h eq "BUSCO" ){ $type{$h}->{$list{$a}->{score}}->{$list{$a}->{$h}}=1;  $tot{$h}->{$list{$a}->{$h}} =1; } 
					else {$type{$h}->{$list{$a}->{score}}++; $tot{$h}++; }
				}
			} else { $out1 .= "NA\t"; }
		}
		print {$fo} "$out1\n";	
	}
}
close($fo);
my %all;
my @ids = ('score', 'tot_rsem_tpm', 'BUSCO');

open(my $fprot, "<", "$res_dir/$genome_id" .".prot.fasta");
open(my $fout, ">", "$res_dir/$out_id" .".prot.fasta");
while(my $v = <$fprot>){

}

foreach my $a (keys(%uni)) {
	my $out1 = $a; foreach my $b (@ids) {
		if ($b ne "BUSCO") {
			if ($a eq "GOOD" || $a eq "BEST" || $a eq "LONGORF" || $a eq "BEST2") { $all{$b} += ($type{$b}->{$a}); }
			if ($tot{$b} > 0) { $out1 .= "\t" . int(10000 *($type{$b}->{$a}/$tot{$b}))/100; } else { $out1 .= "\tNA"; } 
		}
		else {
		 	if ($a eq "GOOD" || $a eq "BEST" || $a eq "LONGORF" || $a eq "BEST2") { 
				foreach my $c (keys(%{$type{$b}->{$a}})) { $all{$b}->{$c} = 1; } 
			}
                        if ($tot{$b}) { $out1 .= "\t" . int(10000 *scalar(keys(%{$type{$b}->{$a}}))/scalar(keys(%{$tot{$b}})))/100; } else { $out1 .= "\tNA"; }
		}
	}
	print STDERR $out1, "\n";
}
 my $out1 = "\nPASS"; foreach my $b (@ids) {
		 if( $tot{$b} > 0) {
               		if ($b ne "BUSCO") { $out1 .= "\t" . int(10000 *($all{$b}/$tot{$b}))/100; } 
	 		else { $out1 .= "\t" . int(10000 *scalar(keys(%{$all{$b}}))/scalar(keys(%{$tot{$b}})))/100;   } }
		 else { $out1 .= "\tNA"; }
	}
        print STDERR $out1, "\n";
sub read_all_rsem($$$)
{
	my $matrix = $_[0];
    my $rsem_file = $_[1];
    my $id = $_[2];
	open( my $fo, "<", $rsem_file) or die("No $rsem_file");
	my $head;
	while (my $v = <$fo> )
    {
		$v =~ /([^\n\r]+)/; my @n = split "\t", $1;
                if (!$head)
		{
			for (my $i = 1; $i < scalar(@n); $i++)
			{
				$head->[$i] = $n[$i];
			}
		}
		else
		{
			if ($n[0] =~ /(\S+)/)
			{
				my $gene_id = $1;
				my $allele = $2;
                		#print STDERR "$gene_id $allele ", $n[0], " $allele ", $matrix->{$gene_id}, "\n";
               			if ($matrix->{$gene_id})
               			{
					my %tmp;
					for (my $col = 1; $col < scalar(@n); $col++)
					{
                    				$matrix->{$gene_id}->{"tot_".$id} += $n[$col];
                        			if (!$tmp{substr($head->[$col], 0, 3)}) { $tmp{substr($head->[$col], 0, 3)}->[0] = $n[$col]; }
						else { push @{$tmp{substr($head->[$col], 0, 3)}}, $n[$col] };
						if ($n[$col] > $matrix->{$gene_id}->{"best_" . $id})
                        			{
                        				$matrix->{$gene_id}->{"best_".$id} = $n[$col];
                        				$matrix->{$gene_id}->{"best_".$id."_id"} = $head->[$col];
						}
                    			}
					my @meds; my $best_med; my $best_lib;
					foreach my $t (keys(%tmp)) {
						my @s = sort {$a <=> $b} @{$tmp{$t}};
						if (scalar(@s) > 2) {
							my $m = -1;
							if (scalar(@s) % 2 == 1) { $m = $s[int(scalar(@s)/2)]; }
							if (scalar(@s) % 2 == 0) {  $m = ($s[int(scalar(@s)/2)-11] + $s[int(scalar(@s)/2)]) /2; }
							if ($m > $best_med) { $best_med = $m; $best_lib = $t; }
						}
					}
		    			$matrix->{$gene_id}->{"avg_" . $id} = $matrix->{$gene_id}->{"tot_".$id} / ( scalar(@n) -1);
					$matrix->{$gene_id}->{"best_median_". $id } = $best_med;
					$matrix->{$gene_id}->{"best_median_lib_". $id } = $best_lib;
                		}
			}

		}
	}
	return($matrix);
}
sub max($$)
{
	if ($_[0] > $_[1]) { return $_[0]; } return ($_[1]); }
sub min($$)
{
        if ($_[0] < $_[1]) { return $_[0]; } return ($_[1]); }


