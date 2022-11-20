#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;

my ($out_dir, $busco5, $qsub, $trinity, $genome_id, $eval, $skip_bowtie, $b_eval, $start, $test, $help, $all_rsem, $blast, $res_dir,$read_dir, $read_end, $mem, $cpu);
$genome_id = "GEN";;
$eval= 1e-5;
$mem = 40;
$cpu = 1;
GetOptions("p|cpu:s"=>\$cpu, "b|busco5"=>\$busco5, "q|qsub"=>\$qsub, "r|result:s"=>\$res_dir,"x|skip"=>\$skip_bowtie, "m|mem:s"=>\$mem,"e|end:s"=>\$read_end, "g|genome:s"=>\$genome_id,  "t|read_dir:s"=>\$read_dir, "o|out_dir:s"=>\$out_dir, "h|?|help"=>\$help);

if ($help || !$read_dir)
{
    print STDERR "Runs Annoation and Translation on a Trimmed Sample\n\n";
    print STDERR "-----------------USAGE--------------\n";
    print STDERR " -t read file directory\n";    
    print STDERR " -e read file ending. Default is \n";
    print STDERR " -r results directory\n";
    print STDERR " -o output directory. Default is the current directory\n";
    print STDERR " -g genome id. Default = GEN\n";
    print STDERR " -m memory for the run. Default = 40g\n";
    print STDERR " -q runs PBS submissions\n";
    print STDERR " -x skip bowtie2 runs\n";
    print STDERR " -b busco is version 5\n";
    print STDERR " -p number of cpus to use per blast submit. Default is 1\n";
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

my $src = Cwd::abs_path(__FILE__);
#Make Databases
my $database_id;
if (-e $res_dir. "/" . $genome_id .".all_nuc.fasta")
{
	open(my $f_o, ">", $res_dir. "/" . $genome_id .".all_nuclear.fasta");
	open(my $f_f, "<", $res_dir. "/" . $genome_id .".all_nuc.fasta");
	while(my $v = <$f_f>){
		while($v =~ />(\S+)/){
			my $id = $1; my $tmp; $v = <$f_f>;
			while ($v !~ />/ && !eof($f_f)){
				$v =~ /([^\n\r]+)/; $tmp .= $1; $v = <$f_f>;
			}
			if (eof($f_f)) {  $v =~ /([^\n\r]+)/; $tmp .= $1; }
			if (length($tmp) > 0 && $id ne $tmp) { print {$f_o} ">$id\n$tmp\n"; }
		}	
	}
	$database_id =  $genome_id . "_ALL";
	my $str = $config_hash->{bowtie} . "\n". $config_hash->{rsem} ."
	cd $res_dir
	bowtie-build " . $res_dir. "/" . $genome_id .".all_nuclear.fasta ".$res_dir ."/" . $database_id ."
	rsem-prepare-reference " . $res_dir. "/" . $genome_id .".all_nuclear.fasta " . $genome_id . "_ALL";
	if (! -e   $res_dir. "/" . $database_id. ".ti")
	{ print STDERR "Making RSEM Database at " .$res_dir. "/" .$database_id ."\n"; `$str`; }
	$str = "cd $res_dir\n" . $config_hash->{bowtie2}."\nbowtie2-build " . $res_dir. "/" . $genome_id .".all_nuc.fasta " . $database_id;
	if (!-e $res_dir. "/" .$database_id . ".1.bt2" && !$skip_bowtie)
	{ print STDERR "Making Bowtie2 Database at " .$res_dir. "/" .$database_id."\n"; `$str`; }
}
else
{
	if (-e $res_dir. "/" . $genome_id .".all_nuc.fasta")
	{
		my $str = $config_hash->{bowtie} . "\n" . $config_hash->{bowtie2} . "\n" . $config_hash->{rsem} ."
		cd $res_dir
		bowtie2-build " . $res_dir. "/" . $genome_id .".all_nuc.fasta ". $genome_id . "
		bowtie-build " . $res_dir. "/" . $genome_id .".all_nuc.fasta ". $genome_id . "
		rsem-prepare-reference " . $res_dir. "/" . $genome_id .".all_nuc.fasta ". $genome_id;
		`$str`;
	}
	else
	{
		die("Cannot find all_nuc data")
	}
}
if (!-e $res_dir . "/rsem")
{
	`mkdir $res_dir/rsem`;
}
if (!-e $res_dir . "/bowtie")
{
	`mkdir $res_dir/bowtie`;
}
my $add;
if ($qsub) { $add = " -q"; }
my $str = $config_hash->{perl} . "\ncd $res_dir/\nperl $src/run_multiple_readexp.pl -d $read_dir -o $out_dir/rsem -p -i $res_dir/$database_id -t rsem -e $read_end -m $mem" . $add;
`$str`;
if (!$skip_bowtie)
{
	$str = $config_hash->{perl} . "\ncd $res_dir/\nperl $src/run_multiple_readexp.pl -d $read_dir -o $out_dir/bowtie -p -i $res_dir/$database_id -t bowtie2 -e $read_end -m $mem" . $add;
	`$str`;
}

my $out = "#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=CPU
#SBATCH --cpus-per-task=CPU
#SBATCH --output=GEN_v_Swissprot.out
#SBATCH --mem=1g
#SBATCH --job-name=JOBNAME
#SBATCH -p intel

";
my $cmd_run = "sbatch";

if ($qsub){
        $cmd_run  = "qsub";
        $out = "#!/bin/bash -l

#PBS -l nodes=1:ppn=CPU
#PBS -l mem=1gb
#PBS -o GEN_v_Swissprot.out
#PBS -e GEN_v_Swissprot.err
#PBS -N JOBNAME
";
}

$out .=
$config_hash->{blast} . "\n" . 
$config_hash->{uniprot} . "

cd ~/
blastp -query ".$res_dir. "/" . $genome_id .".prot.fasta -db \$UNIPROT_DB/uniprot_sprot.fasta -out " . $out_dir. "/" . $genome_id . "_v_SwissProt.blastx.txt -max_target_seqs 1 -num_threads $cpu -evalue 1e-5 -outfmt 6";

$out =~ s/JOBNAME/swissprot_blast/g; $out =~ s/GEN/$genome_id/g; $out =~ s/CPU/$cpu/g;
open(my $fo, ">", $out_dir ."/". $genome_id . "_v_swissprot.sh"); print {$fo} $out; close($fo);
$str = "chmod a+x " . $out_dir ."/". $genome_id . "_v_swissprot.sh\ncd $out_dir\n$cmd_run  ". $out_dir ."/". $genome_id . "_v_swissprot.sh";
if (!-e  $out_dir. "/" . $genome_id . "_v_SwissProt.blastx.txt")
{
	`$str`;
}

$out = "#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=CPU
#SBATCH --output=run_GEN_Trin_busco2_2.txt
#SBATCH --mem=20g
#SBATCH --time=24:00:00
#SBATCH --job-name=run_busco-GEN
#SBATCH -p intel

";
if ($qsub){
        $out = "#!/bin/bash -l
#PBS -l nodes=1:ppn=CPU
#PBS -l mem=20gb
#PBS -o GEN_v_busco.out
#PBS -e GEN_v_busco.err
#PBS -N JOBNAME
";
}


$out .= $config_hash->{busco} ."

cd $out_dir
";
if (!$busco5){
	$out .= "run_BUSCO.py -i ". $res_dir. "/" . $genome_id .".prot.fasta -o $genome_id  -l $blast/arthropoda_odb9 -m prot -c 1 -f -sp tick";
}
else{
	$out .= "busco -i ". $res_dir. "/" . $genome_id .".prot.fasta -o $genome_id  -l arthropoda_odb10 -m prot -c 1 -f";
}
$out =~ s/JOBNAME/busco/g; $out =~ s/GEN/$genome_id/g; $out =~ s/CPU/$cpu/g;
if (!-e  $out_dir . "/run_" . $genome_id)
{
	open(my $fo, ">", $out_dir ."/". $genome_id . "_v_BUSCO.sh"); print {$fo} $out; close($fo);
	$str = "chmod a+x " . $out_dir ."/". $genome_id . "_v_BUSCO.sh\ncd $out_dir\n$cmd_run  ". $out_dir ."/". $genome_id . "_v_BUSCO.sh";
	`$str`;
}
$out = "#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=CPU
#SBATCH --cpus-per-task=1
#SBATCH --output=GEN_v_HMM.out
#SBATCH --mem=1g
#SBATCH --job-name=hmm
#SBATCH -p intel
";

if ($qsub){
        $out = "#!/bin/bash -l
#PBS -l nodes=1:ppn=CPU
#PBS -l mem=1gb
#PBS -o GEN_v_HMM.out
#PBS -e GEN_v_HMM.err
#PBS -N JOBNAME
";
}


$out .= $config_hash->{hmmer} . "\n" .
$config_hash->{pfam} . "

cd $out_dir
hmmscan --noali --tblout $genome_id.hmm.txt -E 1e-5 \$PFAM_DB/Pfam-A.hmm $res_dir"."/"."$genome_id.prot.fasta";
$out =~ s/JOBNAME/busco/g; $out =~ s/GEN/$genome_id/g; $out =~ s/CPU/$cpu/g;

if (!-e  $out_dir . "/$genome_id.hmm.txt")
{
        open(my $fo, ">", $out_dir ."/". $genome_id . "_v_PFAM.sh"); print {$fo} $out; close($fo);
        $str = "chmod a+x " . $out_dir ."/". $genome_id . "_v_PFAM.sh\ncd $out_dir\n$cmd_run  ". $out_dir ."/". $genome_id . "_v_PFAM.sh";
        `$str`;
}

