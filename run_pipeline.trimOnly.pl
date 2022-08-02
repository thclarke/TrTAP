#!/usr/bin/perl
use strict;
use Bio::TreeIO;
use Getopt::Long;
use Bio::TreeIO::NewickParser;
use Cwd;
use Bio::Seq;
use File::Basename;

my ($out_dir, $trinity, $chimera, $genome_id, $eval, $b_eval, $start, $test, $help, $all_rsem, $blast, $res_dir,$closest, $full_chim, $new_genome_id);
$genome_id = "GEN";
$eval= 1e-5;
GetOptions("x|chimera:s"=>\$chimera, "r|result:s"=>\$res_dir, "s|species:s"=>\$b_eval, "e|eval:s"=>\$eval, "a|all_rsem:s"=>\$all_rsem,"g|genome:s"=>\$genome_id,  "t|trinity:s"=>\$trinity, "o|out_dir:s"=>\$out_dir, "b|blast_dir:s"=>\$blast, "c|closest:s"=>\$closest, "z|chimera2:s"=>\$full_chim, "n|new_genome_id:s"=>\$new_genome_id, "h|?|help"=>\$help);


if ($help || !$trinity)
{
    print STDERR "Runs Annoation and Translation on a Trinity Sample\n\n";
    print STDERR "-----------------USAGE--------------\n";
    print STDERR " -t Trinity Assembly (Required)\n";
    print STDERR " -b blast directory\n";
    print STDERR " -r results directory\n";
    print STDERR " -c closest genome. Required";
    print STDERR " -o output directory. Default is the current directory\n";
    print STDERR " -g genome id. Default = GEN\n";
    print STDERR " -e evalue cutoff. Default is 1e-5\n";
    print STDERR " -s species evalue cutoff. Default is the e-value\n";
    print STDERR " -x chimera file. Default is expected _CHIMERA.txt\n";
    print STDERR " -z full chimeric file. No default\n";
	exit();
}

my $config_file = dirname($0) . "/trtap.ini";
if (!-e $config_file){
	print STDERR "Cannot locate config file at $config_file... Exiting...\n\n"; 
	quit(-1);
}

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


if (!$b_eval)
{
	$b_eval = $eval;
}
if (!$new_genome_id){
	$new_genome_id = $genome_id;
}
my %map;
my %seq;
my $matrix;
my $best;
my $chim;
my $chim_num = 1;
my $ls;
my $longest_seq;
my $all_alleles;
if (!$chimera)
{
	$ls = `ls $res_dir/*_CHIMERA.out`;
}
else
{
	$ls = $chimera;
}
while ($ls =~ /([^\n\r]+)/g)
{
	my $file = $1;
	print STDERR "Reading in Chimeric sequences from $file...\n";
	open(my $fc, "<", $file);
	while(my $v = <$fc>)
	{
		while ($v =~ /([^\n\r]+)/g)
		{
			$chim->{$1} = $chim_num;
		}
	}
	$chim_num++;
	print STDERR "Added in ", scalar(keys(%$chim))," contigs as chimeric...\n";
}
my $orig_chim;
if ($full_chim && -e $full_chim)
{
    open(my $fc, "<", $full_chim);
      while(my $v = <$fc>)
      {
                while ($v =~ /([^\n\r]+)/g)
                {
		    $orig_chim->{$1} = $chim_num;
                }
      }

	print STDERR "FOund ", scalar(keys(%$orig_chim))," contigs as chimeric using other file...\n";
}

my $ls1 = "ls $res_dir/$genome_id" . "_v_*RNA.out";
my $ls = `$ls1`;
while ($ls =~ /([^\n\r]+)/g)
{
  	my $file = $1;
        print STDERR "Reading in RNA sequences from $file...\n";
        open(my $fc, "<", $file);
        while(my $v = <$fc>)
        {
                while ($v =~ /([^\n\r]+)/g)
                {
			my $v = $1; $v =~ /(\S+)/;
                        $chim->{$1} = $chim_num;
                }
        }
 	$chim_num++;
}

if ($trinity && -e $trinity)
{
	print STDERR "Reading in $trinity\n";
	open(my $fo, "<", $trinity);
	while(my $v = <$fo>)
	{
		while ($v =~ />(\S+)/ && !eof($fo))
		{
			my $id = $1; my $tmp; $v= <$fo>;
			while($v !~ />/ && !eof($fo))
			{
				$v =~ /([^\n\r]+)/; $tmp .= $1; $v = <$fo>;
			}
			if (eof($fo))
			{
				$v =~ /([^\n\r]+)/; $tmp .= $1;
			}	
			if ($map{$id})
			{
				print STDERR "Warning: $id already seen...\n";
			}
			$map{$id}->{len} = length($tmp);
			$seq{$id} = $tmp;
		}
	}
	print STDERR "Finished reading.. fasta\n";
	foreach my $a (keys(%seq))
	{
		#print "$a\n";
		my $gene_id; my $allele;
		if ($a =~ /(\S+)_i(\d+)\Z/){
			$gene_id = $1;
			$allele = $2;
		}
		else{
			$gene_id = $a;
			$allele = 1;
		}
		$all_alleles->{$gene_id} .= "$allele;";
		if (!$matrix->{$gene_id})
		{
			$matrix->{$gene_id}->{score} = "NO_HIT";
		}
		if (length($seq{$a}) > $matrix->{$gene_id}->{long_len})
		{
			$matrix->{$gene_id}->{long_len} = length($seq{$a});
			$matrix->{$gene_id}->{long_allele} = $allele;
		}
	}
	print STDERR scalar(keys(%$matrix)), " genes identified...\n";
}
else
{
	die("Cannot find trinity file $trinity\n");
}
my $rank;
if (-e $blast)
{
	$ENV{'BLASTDB'} =  $ENV{'BLASTDB'} . ":". $blast;
	#export BLASTDB=$BLASTDB:/bigdata/hayashilab/shared/ayoublab/toby/database/
	my $ls = `ls $blast/*.?hr`;
	while($ls =~ /([^\n\r]+)/g)
	{
		print STDERR "reading $1 blast file..\n";
		my $file = $1;
		if ($file =~ /\A(\S+)\.([pn])hr/)
		{
			my $db_id=$1;
			my $type = $2;
			my $cnt  =0;
			my $dbtype = "prot";
			if ($type eq "n") { $dbtype = "nucl"; } 
			$db_id =~ s/$blast//g;
			while ($db_id =~ /\A\//) { $db_id =~ s/\A\///; }  
			my $out_str = $config_hash->{blast} . "; blastdbcmd -db $db_id -entry all -dbtype $dbtype -outfmt %t--%l";
			my $out = `$out_str`;
			#print "$out_str\n";
			while ($out =~ /([^\n\r]+)/g)
			{
				my $v = $1; $v =~ /(.*)--(.*)/; my @n; $n[0] =$1; $n[1] = $2;  $n[0] =~ /\A(\S+)/; my $id = $1;
				if ($map{$id})
				{
					print STDERR "Warning: $id already seen...\n";
				}
				$cnt++;
				$map{$id}->{len} = $n[1];
			}

			print STDERR "Read in $cnt lines..\n";
		}
	}
	if (-e "$blast/ranking.txt")
	{
		open(my $fo, "<", "$blast/ranking.txt");
		while (my $v = <$fo>)
		{
			$v =~ /([^\n\r]+)/; $v = $1;
			if ($v =~ /\(/)
			{

				my $r = int(rand(1000));
				open(my $fout, ">", "$r");
				print {$fout} $v, "\n";
				close($fout);
				my $treeio = Bio::TreeIO->new(-file => $r, -format=> "newick");
				my $tree = $treeio->next_tree;
				print $tree;
				my $root = $tree->get_root_node;
				my @taxa = $tree->get_leaf_nodes;
				my %ranks;
				`rm $r`;
				if ($closest)
				{
					my @nodes = $tree->find_node(-id => $closest);
					if (scalar(@nodes) < 1)
					{
						print STDERR "Cannot find the node... randomly ranking the nodes...";
						foreach my $a (@taxa)
						{
							$ranks{$a->id} = 1;
						}
					}
					else
					{
						foreach my $a (@taxa)
						{
							$ranks{$a->id} = $tree->distance(-nodes => [$nodes[0], $a] );
						}
					}
				}
				else
				{
					foreach my $a (@taxa)
					{
						$ranks{$a->id} = 1;
					}
				}
				my @sort = sort { $ranks{$a} <=> $ranks{$b} } keys(%ranks);
				push @{$rank}, @sort;

			}
			else
			{
				#print $1, "\n";
				$v =~ /([^\n\r+]+)/; $rank->[scalar(@$rank)] = $1;
			}
		}
	}
	else
	{
		print STDERR "Cannot rank the genomes... quitting."
	}
}

print STDERR "Finished the ranking... $res_dir\n";
if (-e $res_dir)
{
	print STDERR "A\n";
	my $ls1 = "ls $res_dir/$genome_id" . "_v_*.*blastx.out";
	my $ls = `$ls1`;
	while ($ls =~ /([^\/]+)_v_([^\.]+)\.([.]*)blastx.out/g)
	{
		my $genome_id = $1; my $db_id = $2;
		if (!-e $res_dir . "/" . $genome_id ."_" . $db_id .".finished")
		{
			print STDERR "run of $genome_id against $db_id is not finished....\n";
		}
	}
	print STDERR "Finished loading..\n";
	foreach my $a (@$rank)
	{
		if (length($a) > 0)
		{
		        my $ls1 = "ls $res_dir/$genome_id" . "_v_$a.*blastx.out";
			my $ls = `$ls1`;
			if ($ls !~  q/No such file or directory/)
			{
				$ls =~ /([^\n\r]+)/;
				my $file = $1;
				print STDERR "Reading in $file...\n";
				if ($file =~ /(\w+)_v_(\w+)\.(\S*)blastx.out\Z/)
				{
					my $ret = read_blast_file($matrix, $file, $a, \%map, $best, $eval, $chim);
					$matrix = $ret->[0];
					$best = $ret->[1];
				}
			}
		}
	}
	print STDERR "...Finished Loading ", scalar(keys(%$matrix)), " genes.  Starting the analysis...\n";
	my $cnt2;
	my $cnt1;
	my $max_len;
	foreach my $a (keys(%$matrix))
	{
		if (!$matrix->{$a}->{id} && $matrix->{$a}->{score} eq "NO_HIT")
		{
			$cnt2++;
			if ($cnt2 % 10000 == 0)
			{
				print STDERR $cnt2, "\n";
			}
			#print STDERR "Running $a...\n";
			my $orf;
			my @n_1 = split ";", $all_alleles->{$a};
			
			my $bst_c; 	
			foreach my $c (@n_1){
				if ($seq{$a . "_i" . $c} && !$chim->{$a . "_i" . $c}){
					my $tmp = find_longest_orf($seq{$a . "_i" . $c}, 0, 0);
					if (!$orf || length($tmp->{s}) > $orf->{s}){
						$orf = $tmp;
						$bst_c = $c;
					}
					if ($c > 1000){
						die($a.  "_i$c");
					}
				}
				else{
					if ($seq{$a} && !$chim->{$a}){
                                       		 my $tmp = find_longest_orf($seq{$a}, 0, 0);
                                        	if (!$orf || length($tmp->{s}) > $orf->{s}){
                                                	$orf = $tmp;
                                                	$bst_c = $c;
                                     		}
					}
				}
				
			}
			if ($orf)
			{
				if (length($orf->{s}) > $max_len) { $max_len = length($orf->{s}); }
				if (length($orf->{s})  >= 50)
				{
					$matrix->{$a}->{score} = "LONGORF";
				}
				$matrix->{$a}->{id} = $bst_c;
				$matrix->{$a}->{match_genome} = "NA";
				$matrix->{$a}->{match_gene} = "NA";
				$matrix->{$a}->{match_gene_cov} = "NA";
				$matrix->{$a}->{e_value} = "NA";
				$matrix->{$a}->{start} = $orf->{f};
			}
		}
		else
		{
			 $cnt1++;
                        if ($cnt1 % 10000 == 0)
                        {
                                print STDERR $cnt1, "--hits\n";
                        }
		}
	}
	print STDERR "Getting open reading frames..\n";
	my $cnt_m;
	foreach my $a (keys(%$matrix))
	{
		my $g = $a; if ($a =~ /\A2-(\S+)/) { $g = $1; }
		my $gn = $g . "_i" . $matrix->{$a}->{id}; if (!$seq{$gn}) { $gn = $g; }
		if (!($matrix->{$g}->{score} eq "CHIMERA" ||  $matrix->{$g}->{score} eq "RNA"))
		{
			if ($matrix->{$g}->{id} && $seq{$gn})
			{
				$cnt_m++;
				my $tmp_orf = find_longest_orf($seq{$gn}, $matrix->{$a}->{start}, $matrix->{$a}->{dir});
				$matrix->{$a}->{seq} = $tmp_orf->{s};
				$matrix->{$a}->{seq_len} = length($tmp_orf->{s});
				$matrix->{$a}->{nseq} = $tmp_orf->{n};
				$matrix->{$a}->{seq_type} = $tmp_orf->{t};
				$matrix->{$a}->{orig_len} = $map{$gn}->{len};
			}
			else
			{
			 	if (!$matrix->{$a}->{id} && $matrix->{$a}->{score} eq "NO_HIT")
                		{	
                        		$cnt2++;
                        		if ($cnt2 % 10000 == 0)
                        		{
                                		print STDERR $cnt2, "\n";
                        		}
                        		#print STDERR "Running $a...\n";
                        		my $orf;
                        		my @n_1 = split ";", $all_alleles->{$a};

                        		my $bst_c;
                        		foreach my $c (@n_1)
                            		{
                                		if (!$chim->{$gn})
                                		{
                               				 my $tmp = find_longest_orf($seq{$gn}, 0, 0);
                              				  if (!$orf || length($tmp->{s}) > $orf->{s})
                               				 {
                                        			$orf = $tmp;
                                        			$bst_c = $c;
                               			 	}
                              

                        			}
					}
                        		if ($orf)
                        		{
                                		if (length($orf->{s}) > $max_len) { $max_len = length($orf->{s}); }
                                		if (length($orf->{s})  >= 50)
                                		{
                                        		$matrix->{$a}->{score} = "LONGORF";
                                		}
                                		$matrix->{$a}->{id} = $bst_c;
                                		$matrix->{$a}->{match_genome} = "NA";
                                		$matrix->{$a}->{match_gene} = "NA";
                                		$matrix->{$a}->{match_gene_cov} = "NA";
                                		$matrix->{$a}->{e_value} = "NA";
                                		$matrix->{$a}->{start} = $orf->{f};
                        		}
                		}
			}
		}
	}
	print STDERR "Found $cnt_m genes lengths\n";

	print STDERR "Sorting all blast hits..\n";
	$matrix = check_double($matrix, $best);
	print STDERR "Removing Bacterial hits\n";
	my $ls = `ls $res_dir/NR_v_*.txt`;
        my $list = "";
        my %NRmap;
	my $nr_hit = 0;
        while ($ls =~ /NR_v_(\S+).txt/g)
        {
                if ($1 ne "TRINITY")
                {
			$nr_hit++;
                        open(my $fo, "<", "$res_dir/NR_v_$1.txt");
                        while (my $v = <$fo>)
                        {
                                $v =~ /([^\n\r]+)/; my @n = split "\t", $1;
                                if (!$b_eval || $n[10] < $b_eval)
				{
					$NRmap{ids}->{$n[0]} = $n[1];
                                	$NRmap{good}->{$n[1]} = 1;
                        	}
			}
                }
        }
	my @k =  keys(%{$NRmap{good}});
	my $cnt;
	#Getting the orfs of 2nd bests-
	foreach my $a (keys(%$matrix))
	{
		if ($a =~ /\A2-(\S+)/) { 
			my $g = $1; 
		   if ($matrix->{$a}->{id} && $seq{$g . "_i" . $matrix->{$a}->{id}})
                        {
                                $cnt_m++;
				my $gn = $g . "_i" . $matrix->{$a}->{id};
                                if (!$seq{$gn}) { $gn = $a; }
                                my $tmp_orf = find_longest_orf($seq{$gn}, $matrix->{$a}->{start}, $matrix->{$a}->{dir});
                                $matrix->{$a}->{seq} = $tmp_orf->{s};
                                $matrix->{$a}->{seq_len} = length($tmp_orf->{s});
                                $matrix->{$a}->{nseq} = $tmp_orf->{n};
                                $matrix->{$a}->{seq_type} = $tmp_orf->{t};
                                $matrix->{$a}->{orig_len} = $map{$gn}->{len};
                        }
		}
	}
	
	print STDERR "Printing output files..\n";
	print STDERR  $out_dir . "/" . $new_genome_id .".prot.fasta\n";
	open(my $ff, ">", $out_dir . "/" . $new_genome_id .".prot.fasta");
	open(my $fc, ">", $out_dir . "/" . $new_genome_id. ".cds.fasta");
	open(my $fn, ">", $out_dir . "/" . $new_genome_id. ".nuc.fasta");

	open(my $fa, ">", $out_dir . "/" . $new_genome_id. ".all_nuc.fasta");

	open(my $fb, ">", $out_dir . "/" . $new_genome_id. ".bacteria.fasta");
	    open(my $fm, ">", $out_dir . "/" . $new_genome_id . ".gene_info.txt");
	print {$fm} "GeneID\tPipelineType\tBestAllele\tBestGenomeDB\tBestMatch\tMatchLength\tGeneCoverage\tEval\tStartPos\tOrigLength\tFrame\tDirection\tPossible_Chimeric\tSecondary_Allele\tBit_Score\tLength of ongest Allele\tLongest Allele\tProtein_Length\n";
	print STDERR "Running through ", scalar(keys(%$matrix)), " genes...\n";
	    foreach my $a (keys(%$matrix))
	    {
		my $g = $a;
 		my $gn = $a . "_i" . $matrix->{$a}->{long_allele};
                if (!$seq{$gn}) { $gn = $a; }
		if ($a =~ /\A2-(\S+)/) { $g = $1; }#die("$a $g " . $matrix->{$a}->{id}); }
		    if (!$matrix->{$a}->{id})
		    {
			if ($matrix->{$a}->{long_allele})
			{
		    		#my $gn = $a . "_i" . $matrix->{$a}->{long_allele};
				#if (!$seq{$gn}) { $gn = $a; }
				my $orf = find_longest_orf($seq{$gn}, 0,0);
                              	print {$fm} print_no_row($matrix, $a), "\t", length($orf->{s}), "\n";
				print {$fa} ">$a\n", $seq{$gn}, "\n";
			}
			else
			{
				print {$fm} print_no_row($matrix, $a), "\tNA", "\n";
				print STDERR "Cannot find sequence for $a\n";
			}
		    }
		    elsif ($seq{$g . "_i" . $matrix->{$a}->{id}})
		    {
			if ($matrix->{$a}->{species} == "BACTERIA")
			{
				print {$fb} ">$a\n",  $matrix->{$a}->{seq}, "\n";
			}
			if($matrix->{$a}->{score} eq "BEST" || $matrix->{$a}->{score} eq "NOT_BEST" || $matrix->{$a}->{score} eq "LONGORF" || $matrix->{$a}->{score} eq "BEST2")
			{
				if ($orig_chim->{$g . "_i" . $matrix->{$a}->{id}})
				{
				    $matrix->{$a}->{posschi} = 1;
				}
				print {$fm} print_row($matrix, $a), "\t", $matrix->{$a}->{seq_len}, "\n";
				if (!$matrix->{$a}->{species} || $matrix->{$a}->{species} == "EUKARYOTA")
				{
					print {$ff} ">$a\n", $matrix->{$a}->{seq}, "\n";
					print {$fc} ">$a\n", $matrix->{$a}->{nseq}, "\n";
					print {$fa} ">$a\n", $seq{$gn}, "\n";
					print {$fn} ">$a\n", $seq{$gn}, "\n";
				}
			}
		}
		else
			{
                                my $gn = $g . "_i" . $matrix->{$a}->{id};
                                if (!$seq{$gn}) { $gn = $a; }
	        		if ($matrix->{$a}->{species} == "BACTERIA")
                        	{
                                	print {$fb} ">$a\n",  $matrix->{$a}->{seq}, "\n";
                        	}	
                        	if($matrix->{$a}->{score} eq "BEST" || $matrix->{$a}->{score} eq "NOT_BEST" || $matrix->{$a}->{score}){
                                	if ($orig_chim->{$g . "_i" . $matrix->{$a}->{id}})
                                	{
                                	    $matrix->{$a}->{posschi} = 1;
                                	}
                                	print {$fm} print_row($matrix, $a), "\t", $matrix->{$a}->{seq_len}, "\n";
                                	if (!$matrix->{$a}->{species} || $matrix->{$a}->{species} == "EUKARYOTA")
                                	{
                                        	print {$ff} ">$a\n", $matrix->{$a}->{seq}, "\n";
                                        	print {$fc} ">$a\n", $matrix->{$a}->{nseq}, "\n";
                                        	print {$fa} ">$a\n", $seq{$gn}, "\n";
                                        	print {$fn} ">$a\n", $seq{$gn}, "\n";
                                	}
                        	}
				else{
					my $orf = find_longest_orf($seq{$gn}, $matrix->{$a}->{start}, $matrix->{$a}->{dir});
					if ($orig_chim->{$g . "_i" . $matrix->{$a}->{id}})
                			{
                				$matrix->{$a}->{posschi} = "1";
                			}
					print {$fm} print_no_row($matrix, $a), "\t", length($orf->{s}), "\n";
					print {$fa} ">$a\n", $seq{$g . "_i" . $matrix->{$a}->{id}},"\n";
				}
			}
		

	}

}

sub make_orfs($$)
{
	my $m = $_[0];
	my $s = $_[1];
	my $c = 0; 
	foreach my $a (keys(%$m))
	{
		if (!( $m->{$a}->{score} eq "CHIMERA" ||  $m->{$a}->{score} eq "RNA"))
		{
			if ($m->{$a}->{id})
			{
				$c++;
				my $tmp_orf = find_longest_orf($s->{$a . "_i" . $m->{$a}->{id}}, $m->{$a}->{start}, $m->{$a}->{dir});
				$m->{$a}->{seq} = $tmp_orf->{s};
				$m->{$a}->{seq_len} = length($tmp_orf->{s});
				$m->{$a}->{nseq} = $tmp_orf->{n};
				$m->{$a}->{start_pos} = $tmp_orf->{p};
				
			}
		}
	}
	print STDERR "Found $c genes lengths\n";
	return($m);
}

sub print_row($$)
{
	my $m = $_[0];
	my $id = $_[1];
	my $ret = "";
	$ret .= $id;
	$ret .= "\t" . $m->{$id}->{score};
	#$ret .= "\t" . $m->{$id}->{species};
	$ret .= "\t" . $m->{$id}->{id};
	$ret .= "\t" . $m->{$id}->{match_genome};
	$ret .= "\t" . $m->{$id}->{match_gene};
	$ret .= "\t" . $m->{$id}->{match_len};
	$ret .= "\t" . $m->{$id}->{len_cov};
	$ret .= "\t" . $m->{$id}->{e_value};
	$ret .= "\t" . $m->{$id}->{start_pos};
	$ret .= "\t" . $m->{$id}->{orig_len};
	$ret .= "\t" . $m->{$id}->{start};
	$ret .= "\t" . $m->{$id}->{dir};
	$ret .= "\t" . (0+$m->{$id}->{posschi});
	$ret .= "\t" . (0+$m->{$id}->{sec});
	$ret .= "\t" . (0+$m->{$id}->{bit_score});
 	$ret .= "\t" . $m->{$id}->{long_len};
        $ret .= "\t" .$m->{$id}->{long_allele};
	return($ret);
}

sub print_no_row($$)
{
	my $m = $_[0];
	my $id = $_[1];
	my $ret = "";
	$ret .= $id;
	$ret .= "\t" .$m->{$id}->{score};
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" . (0+$m->{$id}->{posschi});
	$ret .= "\t" ."NA";
        $ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	$ret .= "\t" ."NA";
	return($ret);
}




#Make a function that reads in the matrix and the blast file:
sub read_blast_file($$$$$$$)
{
	my $matrix = $_[0];
	my $blast_file = $_[1];
	my $id = $_[2];
	my $map = $_[3];
	my $best = $_[4];
	my $e_val = $_[5];
	my $chim = $_[6];
	my %hits;
	if (-e $blast_file)
	{
		open(my $fo, "<", $blast_file);
		while( my $v = <$fo> )
		{
			$v =~ /([^\n\r]+)/; my @n = split "\t", $1;
			#if ($n[1] eq "GBL86543.1") { print STDERR "here ", $n[11], " ", $n[0], " ", $n[1]," ", $chim->{$n[0]},  "\n"; }

			#get the gene name
			if ($chim->{$n[0]})
			{
				my $gene_id = $n[0];
                                my $allele = 1;
				if ($n[0] =~ /(\S+)_i(\d+)\Z/){
                        		$gene_id = $1;
                        		$allele = $2;
				}
				#die("$gene_id " . $chim->{$n[0]}. " " . $matrix->{$gene_id}->{score});
                        	if (!$matrix->{$gene_id}->{score} || $matrix->{$gene_id}->{score} eq "NO_HIT")
				{
					if ($chim->{$n[0]} == 1)
					{
						$matrix->{$gene_id}->{score} = "CHIMERA";
					}
					else
					{
						$matrix->{$gene_id}->{score} = "RNA"
					}
				}
			}
			else
			{
				if (!$e_val || $n[10] <= $e_val)
				{
					
					my $gene_id = $n[0]; my $allele = 1;
					if ($n[0] =~ /(\S+)_i(\d+)\Z/){
						$gene_id = $1;
						$allele = $2;
					}
					if ($best->{$id . "--" . $n[1]}->{score} < $n[11])
					{
						if (!$map->{$n[0]}->{len}) {die($n[0]); }
						$best->{$id . "--" . $n[1]}->{score} = $n[11];
						$best->{$id . "--" . $n[1]}->{id} = $gene_id;
						$best->{$id . "--" . $n[1]}->{allele} = $allele;
						$best->{$id . "--" . $n[1]}->{match_gene_cov} = abs($n[9]-$n[8]) / $map->{$n[1]}->{len};						
						$best->{$id . "--" . $n[1]}->{match_self_cov} = abs($n[7]-$n[6]) / $map->{$n[0]}->{len};
						$best->{$id . "--" . $n[1]}->{start} = (($n[6]-1) % 3)+1;
						$best->{$id . "--" . $n[1]}->{dir} = 1;
						$best->{$id . "--" . $n[1]}->{match_len} = $map->{$n[1]}->{len};
						$best->{$id . "--" . $n[1]}->{e_value} = $n[10];
						$best->{$id . "--" . $n[1]}->{bit_score} = $n[11];
						if ($n[7] <  $n[6])
						{
							$best->{$id . "--" . $n[1]}->{start} = (($map->{$n[0]}->{len} - $n[6]) % 3)+1;
							$best->{$id . "--" . $n[1]}->{dir} =   -1;
						}
					}
					if (!$matrix->{$gene_id}->{id} && $map->{$n[0]})
					{
						if (!$map->{$n[1]} || $map->{$n[1]}->{len} == 0 ||  $map->{$n[0]}->{len}==0)
                                		{
                                			my $skip;
				        		warn("From ". $gene_id. " and $id hit " .$n[1] . " is empty " . $map->{$n[1]}->{len});
                                		}
						else
						{
							$hits{$gene_id} = 1;
							$matrix->{$gene_id}->{score} = "HIT";
							$matrix->{$gene_id}->{id} = $allele;
							$matrix->{$gene_id}->{match_genome} = $id;
							$matrix->{$gene_id}->{match_gene} = $n[1];
							$matrix->{$gene_id}->{match_len} =  $map->{$n[1]}->{len};
				#if (!$map->{$n[1]} || $map->{$n[1]}->{len} == 0)
				#{
			#		die("From ". $gene_id. " and $id hit " .$n[1] . " is empty " . $map->{$n[1]}->{len});
		#		}
							$matrix->{$gene_id}->{match_gene_cov} = abs($n[9]-$n[8]) / $map->{$n[1]}->{len};
							$matrix->{$gene_id}->{match_self_cov} = abs($n[7]-$n[6]) / $map->{$n[0]}->{len};
							$matrix->{$gene_id}->{bit_score} = $n[11];
							$matrix->{$gene_id}->{e_value} = $n[10];
							$matrix->{$gene_id}->{start} = 1+(($n[6]-1) % 3);
							$matrix->{$gene_id}->{sec} = 1;
							$matrix->{$gene_id}->{dir} = 1;
							if ($n[7] <  $n[6])
							{
								$matrix->{$gene_id}->{start} = (($map->{$n[0]}->{len} - $n[6]) % 3)+1;
								$matrix->{$gene_id}->{dir} =   -1;
							}
						}
					}
					else
					{
						if ( $matrix->{$gene_id}->{match_genome} eq $id && $matrix->{$gene_id}->{bit_score} < $n[11])
						{
							if ($map->{$n[1]}->{len} > 0)
							{


								$hits{$gene_id} = 1;
								$matrix->{$gene_id}->{id} = $allele;
								$matrix->{$gene_id}->{score} = "HIT";
								$matrix->{$gene_id}->{match_genome} = $id;
								$matrix->{$gene_id}->{match_gene} = $n[1];
								$matrix->{$gene_id}->{match_len} =  $map->{$n[1]}->{len};
								$matrix->{$gene_id}->{match_gene_cov} = abs($n[9]-$n[8]) / $map->{$n[1]}->{len};
								$matrix->{$gene_id}->{match_self_cov} = abs($n[7]-$n[6]) / $map->{$n[0]}->{len};
								$matrix->{$gene_id}->{bit_score} = $n[11];
								$matrix->{$gene_id}->{e_value} = $n[10];
								$matrix->{$gene_id}->{start} = (($n[6]-1) % 3)+1;
								$matrix->{$gene_id}->{sec}++;
								$matrix->{$gene_id}->{dir} = 1;
								if ($n[7] <  $n[6])
								{
									$matrix->{$gene_id}->{start} = (($map->{$n[0]}->{len} - $n[6]) % 3)+1;
									$matrix->{$gene_id}->{dir} =   -1;
								}
							}
							else
							{
								print STDERR $n[1].  " does not have length.. please check..\n";
							}
						}
					}
				}
			}
		}
	}
	print STDERR "Run on $id added ", scalar(keys(%hits)), " genes\n";
	my $ret; $ret->[0] = $matrix; $ret->[1] = $best;

	return $ret;
}



sub check_double($$$)
{
	my $matrix = $_[0];
	my $best = $_[1];
	foreach my $a (keys(%$matrix))
	{
		#Find the best hit per gene
		if ($matrix->{$a}->{match_gene} && $matrix->{$a}->{match_len})
		{
				$matrix->{$a}->{len_cov} = $matrix->{$a}->{seq_len}/$matrix->{$a}->{match_len};
		}
		if ( $matrix->{$a}->{score} eq "HIT")
		{
			if ($best->{$matrix->{$a}->{match_genome} . "--" .$matrix->{$a}->{match_gene}}->{id} eq $a)
				{
					$matrix->{$a}->{score} = "BEST";
				}
				else
				{
					my $b = $best->{$matrix->{$a}->{match_genome} . "--" .$matrix->{$a}->{match_gene}}->{id};
					if($matrix->{$a}->{match_len} > 0 && $matrix->{$a}->{score} ne "CHIMERA" && $matrix->{$a}->{score} ne "RNA")
					{ 	
							$matrix->{$a}->{score} = "NOT_BEST";
					}					
					else
					{
						if ($matrix->{$a}->{match_gene} && $matrix->{$a}->{match_gene} ne "NA")
						{
							#die($matrix->{$a}->{match_gene} . " " .$matrix->{$a}->{match_len});
						}
					}
				}
			
		}
	}
	my $b_cnt;
	foreach my $a (keys(%$best))
	{
		my $id = $best->{$a}->{id};
		my $allele = $best->{$a}->{allele};
		my $gene_cov = $best->{$a}->{match_gene_cov};
		my $self_cov = $best->{$a}->{match_self_cov};
		$a =~ /(\S+)--(\S+)/; my $genome = $1; my $gene_id = $2;
		if ($matrix->{$id}->{id} ne $allele && $matrix->{$id}->{match_genome} eq $genome && $gene_cov > 0.8)
		{
			my $id2 = "2-" . $id;
			$matrix->{$id2}->{id} = $allele;
			$matrix->{$id2}->{match_genome} = $genome;
			$matrix->{$id2}->{match_gene} = $gene_id;
			 $matrix->{$id2}->{match_len} = $best->{$a}->{match_len};
			 $matrix->{$id2}->{len_cov} = $gene_cov;
                         $matrix->{$id2}->{e_value} = $best->{$a}->{e_value};
                        $matrix->{$id2}->{bit_score}= $best->{$a}->{bit_score};
			$matrix->{$id2}->{score} = "BEST2";
			$matrix->{$id2}->{start} = $best->{$a}->{start};
			$matrix->{$id2}->{dir} = $best->{$a}->{dir};
			$b_cnt++;
		}
	}
	print STDERR "Missing 2nd best hits: $b_cnt\n";
	return $matrix;
}

sub translater($$$$)
{
	my $u_dir = $_[2];
	my $u_frame = $_[1];
	my $strn = $_[0];
	my $max = $_[3];
	#else { $strn = substr($strn); }
	$strn = substr($strn, $u_frame);
	my $seq_obj = Bio::Seq->new(-seq => $strn, -alphabet => 'dna'); my $prot_seq = $seq_obj->translate(); my $tr = $prot_seq->seq;
	#my $tr = translate_as_string($strn); 
	$tr .= "#"; my $c = 0;
	while ($tr =~ /([^\*\#]*)([\*\#])/g)
	{
		my $s = $1; my $e = $2;
		if (!$max->{s} || length($s) > length($max->{s}))
		{
			$max->{sl} = length($max->{s}); $max->{s} = $s; $max->{c} = $c; $max->{e} = $e; $max->{f} = $u_frame; $max->{b} = 1;
			$max->{p} = pos($tr) - length($s)-1; $max->{d} = $u_dir; $max->{n} = substr($strn,  3 * $max->{p},  (length($s)) * 3);
    			$max->{t} =  2*($e eq"#") + ($c==1);
    		}
		else
		{
			if (length($s) > $max->{sl})
			{
				$max->{sl} = length($s);
			}
		}
		$c++;
	}
	return $max;
}

#Go through the given frame each of the 6 frames, and return the longest
sub find_longest_orf($$$)
{
	my $str = $_[0];
	my $u_dir = $_[2]; #given direction
	my $u_frame = $_[1]-1; #given frame
	my $max;
	if ($u_dir != 0 && $u_frame > -1) #This is the given finding the longest ORF
	{
		$max =  translater($str, $u_frame, $u_dir, $max);
	}
	else
	{
		for (my $i = 0; $i < 3; $i++)
		{
			#my $tr = translate_as_string(substr($str, $i));
			$max = translater($str, $i, 1, $max);
    			$max = translater($str, $i, -1, $max);
    		}
	}
    	return $max;
}
