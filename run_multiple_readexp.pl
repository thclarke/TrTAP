#!/usr/bin/perl
use Getopt::Long;
use strict;
use warnings;

my ($qsub,$dir, $out, $pattern, $help, $index, $nthreads, $end, $paired, $type, $mem);
$end = ".fq";
$type = "bowtie2";
$mem = 20;
$nthreads = 1;
GetOptions("q|qsub"=>\$qsub, "m|mem:s" => \$mem, "d|dir:s" =>\$dir, "n|nthread"=>\$nthreads,"o|out:s"=>\$out, "t|type:s"=>\$type, "e|end:s"=>\$end, "i|index:s"=>\$index, "paired|p"=>\$paired, "h|?|help"=>\$help);

if ($help || !($dir))
{
    print STDERR "Run bowtie or RSEM on the grid\n";
    print STDERR "-----------------USAGE--------------\n";
    print STDERR " -d directory with reads \n";
    print STDERR " -o output directory\n";
    print STDERR " -p is paired reads\n";
    print STDERR " -i comma separated list of indices\n";
    print STDERR " -t program type to run. Default is bowtie2\n";
    print STDERR " -e end of string\n";
    print STDERR " -n number of threads\n";
    print STDERR " -q runs PBS submissions\n";
    print STDERR " -m memory (in g) to request. Default is 20\n";
	die();
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

my $sh = "#!/bin/bash -l\n\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=NTHREAD\n";
$sh .= "#SBATCH --output=run_rsem_OUT.txt\n\n";
$sh .= "#SBATCH --mem=MEMg\n#SBATCH --time=48:00:00\n#SBATCH --job-name=run_rsem_OUT\n#SBATCH -p intel\n\n";

my $cmd_run = "sbatch";

if ($qsub){
        $cmd_run  = "qsub";
        $sh = "#!/bin/bash -l

#PBS -l nodes=1:ppn=CPU
#PBS -l mem=MEMgb
#PBS -o run_rsem_OUT.out
#PBS -e run_rsem_OUT.err
#PBS -N run_rsem_OUT
";
}


if (lc($type) eq "rsem")
{
	$sh .=  $config_hash->{bowtie} . "\n". $config_hash->{rsem} ."\n";
	$sh .= " cd DIR\nrsem-calculate-expression FQ RSEM_INDEX OUTFILE\n";
}
if (lc($type) eq "bowtie2")
{
	$sh .= $config_hash->{bowtie2} . "\n". $config_hash->{python3} ."\n";
	$sh .= "cd DIR\nbowtie2 -a FQ -x RSEM_INDEX -S OUTFILE.sam\n";
	$sh .= "samtools view -bS OUTFILE.sam > OUTFILE.bam\n";
	$sh .= "rm OUTFILE.sam\n";
	$sh .= "samtools sort -o OUTFILE.sort.bam OUTFILE.bam\n";
	$sh .= "rm OUTFILE.bam\n";
	$sh .= "samtools index OUTFILE.sort.bam\n";
	$sh .= "python3 ~/count_bam_hits_fpkm.py OUTFILE.sort.bam 2 >OUTFILE.bowtie.cnts\n"
}
if (lc($type) eq "count")
{
        $sh .= $config_hash->{bowtie2} . "\n". $config_hash->{python3} ."\n";   
        $sh .= "cd DIR\n";
        $sh .= "python3 ~/count_bam_hits_fpkm.py OUTFILE.sort.bam 2 >OUTFILE.bowtie.cnts\n"
}
my $ls;
if ($paired)
{
	$ls = `ls $dir/*1$end`;
}
else
{
	$ls = `ls $dir/*$end`;
}
my @index = split ",", $index;

while ($ls =~ /([^\n\r]+)/g)
{
	my $file = $1;
	my $fq = $file;
	my $good = 1;
	my @n = split "\/", $file;
	my $out_id = $n[scalar(@n)-1];

	$out_id =~ s/$end//;
	if (!-e $fq)
	{
		$good = 0;
	}
	if ($paired)
	{
		$out_id =~ /(\S)1\Z/;
		my $id1 = $1 . "1" . $end;
		my $id2 = $1 . "2". $end;
		$out_id =~ s/1//;
		if (lc($type) eq "rsem")
		{
			$fq = "--paired-end $file";
		}
		else
		{
			$fq = " -1 $file";
		}
		my $file2 = $file;
		$file2 =~ s/$id1/$id2/;
		if (-e $file && -e $file2)
		{
			if (lc($type) eq "rsem")
			{
				$fq .= " $file2";
			}
			else
			{
				$fq .= " -2 $file2";
			}
		}
		else
		{
			$good = 0;
		}
	}
	if ($good)
	{
		foreach my $a (@index)
		{
		my @n = split "\/", $a;
		my $out_id1 = $out_id . "_" . $n[scalar(@n)-1];
		my $tmp_sh = $sh;
		$tmp_sh =~ s/DIR/$dir/g;

		$tmp_sh =~ s/FQ/$fq/g;
		$tmp_sh =~ s/RSEM_INDEX/$a/g;
		my $outfile = $out . "\/" . $out_id1;
		$tmp_sh =~ s/OUTFILE/$outfile/g;
		$tmp_sh =~ s/OUT/$out_id1/g;
		$tmp_sh =~ s/MEM/$mem/g;
		$tmp_sh =~ s/NTHREAD/$nthreads/g;
		print STDERR "$out_id1\n";
		open(my $fo, ">", "$out\/$out_id1". "_$type.sh");
		print {$fo} $tmp_sh;
		close($fo);
		my $cmd = "chmod a+x ".$out ."\/" .$out_id1 . "_$type.sh\n";
		`$cmd`;
		$cmd = "$cmd_run ".$out ."\/" . $out_id1. "_$type.sh";
		`$cmd`;
		}
	}

}

