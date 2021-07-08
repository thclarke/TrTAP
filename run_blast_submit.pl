#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;
#use Config::File;
use File::Basename;

sub rndStr{ join'', @_[ map{ rand @_ } 1 .. shift ] }

sub get_full_path($){
	my $trinity_file = basename($_[0]);
	my $trinity_dir = dirname($_[0]);
	my $full_dir = `cd "$trinity_dir" && pwd `;
	chomp $full_dir;
	return $full_dir. "/". $trinity_file;
}

sub get_full_dir($){
	my $trinity_dir = $_[0];
        my $full_dir = `cd "$trinity_dir" && pwd `;
        chomp $full_dir;
        return $full_dir. "/";
}

my ($qsub, $chim, $out_dir, $trinity, $random, $skip_nr, $email, $test, $help, $blast, $fq_dir, $num, $test, $rna, $time); 
$email = "toby.h.clarke\@gmail.com";
GetOptions("q|qsub"=>\$qsub, "c|chim:s"=>\$chim, "i|id:s"=>\$random, "k|skip"=>\$skip_nr, "e|email:s" =>\$email, "r|rna"=>\$rna, "m|time:s"=>\$time, "x|test" =>\$test, "f|fq:s"=>\$fq_dir, "t|trinity:s"=>\$trinity, "o|out_dir:s"=>\$out_dir, "b|blast_dir:s"=>\$blast, "n|num:s"=>\$num, "h|?|help"=>\$help);

if ($help || !$trinity)
{
    print STDERR "Runs the Pilon on a reduced Sample\n\n";
    print STDERR "-----------------USAGE--------------\n";
    print STDERR " -t Trinity Assembly (Required)\n";
    print STDERR " -b blast directory (Required)\n";
    print STDERR " -e email address. Default = 0\n";
    print STDERR " -o output directory. Default is the current directory\n";
    print STDERR " -n number to split into. Default = 0\n";
    print STDERR " -x run a test. Does not submit . Default = 0\n";
    print STDERR " -f fastq directory to run paired end RSEM.\n"; 
    print STDERR " -i for the output\n";
    print STDERR " -c database to run the chimera test on\n";
    print STDERR " -q runs PBS submissions\n";
    print STDERR " -k skip nr run\n";
    print STDERR " -r also runs rrna and trna searchs\n";
    print STDERR " -m add time to run\n";
    exit();
}

if (!-e $trinity)
{
	warn "Cannot find $trinity file... quitting";
	quit(-1);
}

$trinity = get_full_path($trinity);
$blast = get_full_dir($blast);
$fq_dir = get_full_dir($fq_dir);
$out_dir = get_full_dir($out_dir);

if ($num)
{
	`perl /rhome/tclarke/split_fasta_file.pl $trinity $num`;
}

my $cmd_run = "sbatch";
my $cmd_multi = "sbatch --array=1-$num";

print $0, "\n";
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

#my $config_hash = Config::File::read_config_file($config_file);

my $ls_nblast =`ls $blast/*nhr`;
my $ls_pblast = `ls $blast/*.phr`;


my $st = "#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=CPU
#SBATCH --output=OUT2.out
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=ALL
#SBATCH --mem=10g
#SBATCH --job-name=JOBNAME
#SBATCH -p intel
";
if ($qsub){
	$cmd_run  = "qsub";
	$cmd_multi = "qsub -J 1-$num";
	$st = "#!/bin/bash -l

#PBS -l nodes=1:ppn=CPU
#PBS -l mem=10gb
#PBS -o OUT2.out
#PBS -e OUT2.err
#PBS -M EMAIL
#PBS -m abe
#PBS -N JOBNAME
";

}
if ($time && !$qsub) { $st .= "#SBATCH --time=". $time . "\n"; }
if ($time && $qsub) { $st .= "#PBS -l walltime=". $time . "\n"; }
my $in = $st . "\n" . $config_hash->{blast} ."

cd DIR
if [  -e RAND_DB.sub ]
then
	rm RAND_DB.sub
fi
touch RAND_DB.NUM.start
BLAST -query TRINITY -db DB1 -out OUT1 -num_threads 1 -evalue 1e-5 -max_target_seqs 1 -outfmt 6
rm RAND_DB.NUM.start
touch RAND_DB.NUM.finished
";

my $in2 = $st . "\n" . 
$config_hash->{rsem} . "\n" .
$config_hash->{bowtie} . "\n
cd DIR

touch RAND_ID.rsemstart
rsem-calculate-expression --paired-end ID1 ID2 DB1 OUT1
rm RAND_ID.rsemstart
if [ -e OUT1.genes.results ]
then

	touch RAND_ID.rsemdone
	rm OUT1.transcript.bam
	touch RAND_ID.rsemclean
else
	touch RAND_ID.error
fi
";

my $in3 = $st . "\n" .
$config_hash->{diamond} . "\n

cd DIR
touch RAND_NR.start
diamond blastx --query TRINITY --db \$DIAMOND_DB/nr --threads 1 --out NR_v_RAND.txt --threads 1 --evalue 1e-5 --max-target-seqs 1  --outfmt 6
rm RAND_NR.start
touch RAND_NR.finished
";

my $in_chim = $st . 
$config_hash->{blast} 
."
".
$config_hash->{python2}
. "

cd DIR
if [  -e RAND_CHIMERA.sub ]
then
	rm RAND_CHIMERA.sub
fi
touch RAND_CHIMERA.NUM.start
BLAST -query TRINITY -db DB1 -out OUT1 -num_threads 1 -evalue 1e-5 -outfmt 6
rm RAND_CHIMERA.NUM.start
if [  -e OUT1 ]
then
	touch RAND_CHIMERA_PY.start
	python2 CURRDIR/src/AnalyzeMultiHits.py -csvIn OUT1 -txtOut OUT2 -overlap 40
	rm RAND_CHIMERA_PY.start
	touch RAND_CHIMERA_PY.finished
fi

touch CHIMERA.NUM.finished
";


my $in_rrna = $st . $config_hash->{blast} 
. "
". $config_hash->{rrna} .
"
cd DIR
if [  -e RAND_RRNA.sub ]
then
        rm RAND_RRNA.sub
fi
touch RAND_RRNA.NUM.start
blastn -max_target_seqs 1 -query TRINITY -db \$RRNADB -out OUT1 -num_threads 1 -evalue 1e-5 -outfmt 6
rm RAND_RRNA.NUM.start
touch RAND_RRNA.NUM.finished
";

my $in_trna = $st . $config_hash->{tRNAscan} 
. "

cd DIR
touch RAND_TRNA.NUM.start
tRNAscan-SE -G -o OUT1 TRINITY
rm RAND_TRNA.NUM.start
touch RAND_TRNA.NUM.finished
";

my @chars = ("A".."Z", "a".."z");my $cpu = 1;
if (!$random) { $random = rndStr 4, @chars; }
if ($chim)
{
	my $type;
	if (-e "$blast/$chim.nhr")
	{
		$type = "n";
		
	}
	if (-e "$blast/$chim.phr")
	{
		$type = "p";
		
	}
	if ($type)
	{
		my $file = "$blast/$chim.". $type . "hr"; $file =~ /\.([np])hr/; my $type = $1; my $b_type;
		$file =~ s/\.([np])hr//; $file =~ /([^\/]+)\Z/; my $id = $1; my $new = $in_chim;
		my $out2 = $out_dir . "/" . $random . "_CHIMERA.out"; $new =~ s/OUT2/$out2/g;
		my $out3 = $out_dir . "/" . $random . "_CHIMERA.out"; $new =~ s/OUT2/$out2/g;
		 $new =~ s/CPU/$cpu/;
		my $file1= $blast . "/" . $id; $new =~ s/DB1/$file1/g; $new =~ s/RAND/$random/g;  my $jobname = "CHIMERA-$random-RUN";
		if ($type eq "p") { $new =~ s/BLAST/blastx/; $b_type = "blastx"; }
		if ($type eq "n") { $new =~ s/BLAST/tblastx/; $b_type = "tblastx";}
		my $abs = dirname($0);
		$new =~ s/CURRDIR/$abs/g;
		$new =~ s/DIR/$out_dir/g;  $new =~ s/EMAIL/$email/; $new =~ s/JOBNAME/$jobname/;
		if (!$num)
		{
			my $out_id = $out_dir ."/" . $random . "_v_CHIMERA.$b_type.out";
			$new =~ s/OUT1/$out_id/g;
			$new =~ s/TRINITY/$trinity/g;
			$new =~ s/NUM\.//g;
		}
		else
		{
			my $out_id = $out_dir ."/" . $random . "_v_CHIMERA.$b_type.\$SLURM_ARRAY_TASK_ID.out";
			$new =~ s/OUT1/$out_id/;
			my $tmp_trinity = $trinity  . ".\$SLURM_ARRAY_TASK_ID";
			$new =~ s/TRINITY/$tmp_trinity/g;
			$new =~ s/NUM/\$SLURM_ARRAY_TASK_ID/g; 
	}
	if (!-e $out_dir . "/" . $random . "_CHIMERA.start" && !-e $out_dir . "/" . $random . "_CHIMERA.finished")
	{  
		open(my $fo, ">", $out_dir . "/" .$random . "_v_CHIMERA.sh");
		print {$fo} $new;
		close($fo);
		my $tmp_sh = $out_dir . "/" .$random . "_v_CHIMERA.sh";
		`chmod a+x $tmp_sh`;
		if (!$test && !-e "$out_dir/". $random ."_CHIMERA.sub") 
		{ 
			my $touch = "touch ". $random ."_CHIMERA.sub";
			`$touch`;

			if (!$num) {`$cmd_run $tmp_sh`; }
			else {	
				`$cmd_multi $tmp_sh`; }
		}
	}
}

}

if ($rna)
{
     	my $new = $in_rrna;
     	my $out2 = $out_dir . "/" . $random . "_RRNA.out"; $new =~ s/OUT2/$out2/g;
	$new =~ s/RAND/$random/g;  my $jobname = "RRNA-$random-RUN";
     	$new =~ s/DIR/$out_dir/g;   $new =~ s/CPU/$cpu/; $new =~ s/EMAIL/$email/; $new =~ s/JOBNAME/$jobname/;
        if (!$num)
	{
      		my $out_id = $out_dir ."/" . $random . "_v_RRNA.out";
                $new =~ s/OUT1/$out_id/g;
                $new =~ s/TRINITY/$trinity/g;
                $new =~ s/NUM\.//g;
        }
        else
       {
        	my $out_id = $out_dir ."/" . $random . "_v_RRNA.\$SLURM_ARRAY_TASK_ID.out";
                $new =~ s/OUT1/$out_id/;
                        my $tmp_trinity = $trinity  . ".\$SLURM_ARRAY_TASK_ID";
                        $new =~ s/TRINITY/$tmp_trinity/g;
                        $new =~ s/NUM/\$SLURM_ARRAY_TASK_ID/g;
        }
        if (!-e $out_dir . "/" . $random . "_RRNA.start" && !-e $out_dir . "/" . $random . "_RRNA.finished")
        {
                open(my $fo, ">", $out_dir . "/" .$random . "_v_RRNA.sh");
                print {$fo} $new;
                close($fo);
                my $tmp_sh = $out_dir . "/" .$random . "_v_RRNA.sh";
                `chmod a+x $tmp_sh`;
                if (!$test && !-e "$out_dir/". $random ."_RRNA.sub")
                {
                        my $touch = "touch ". $random ."_RRNA.sub";
                        `$touch`;
                        if (!$num) {`$cmd_run $tmp_sh`; }
                        else {
                                `$cmd_multi $tmp_sh`; }
                }
        }
 	my $new = $in_trna;
        my $out2 = $out_dir . "/" . $random . "_TRNA.out"; $new =~ s/OUT2/$out2/g;
        $new =~ s/RAND/$random/g;  my $jobname = "TRNA-$random-RUN";
        $new =~ s/DIR/$out_dir/g;  $new =~ s/EMAIL/$email/;  $new =~ s/CPU/$cpu/; $new =~ s/JOBNAME/$jobname/;
        if (!$num)
        {
                my $out_id = $out_dir ."/" . $random . "_v_TRNA.out";
                $new =~ s/OUT1/$out_id/g;
                $new =~ s/TRINITY/$trinity/g;
                $new =~ s/NUM\.//g;
        }
        else
       {
                my $out_id = $out_dir ."/" . $random . "_v_TRNA.\$SLURM_ARRAY_TASK_ID.out";
                $new =~ s/OUT1/$out_id/;
                        my $tmp_trinity = $trinity  . ".\$SLURM_ARRAY_TASK_ID";
                        $new =~ s/TRINITY/$tmp_trinity/g;
                        $new =~ s/NUM/\$SLURM_ARRAY_TASK_ID/g;
        }
        if (!-e $out_dir . "/" . $random . "_TRNA.start" && !-e $out_dir . "/" . $random . "_TRNA.finished")
        {
                open(my $fo, ">", $out_dir . "/" .$random . "_v_TRNA.sh");
                print {$fo} $new;
                close($fo);
                my $tmp_sh = $out_dir . "/" .$random . "_v_TRNA.sh";
                `chmod a+x $tmp_sh`;
                if (!$test && !-e "$out_dir/". $random ."_TRNA.sub")
                {
                        my $touch = "touch ". $random ."_TRNA.sub";
                        `$touch`;
                        if (!$num) {`$cmd_run $tmp_sh`; }
                        else {
                                `$cmd_multi $tmp_sh`; }
                }
        }
}
while ($ls_nblast =~ /([^\n\r]+)/g)
{
	my $file = $1; $file =~ s/\.nhr//; $file =~ /([^\/]+)\Z/; my $id = $1; my $new = $in; my $cpu = 1;
	if ($id eq "Avent") { $cpu = 10; }
	my $out2 = $out_dir . "/" . $random . "_" . $id . ".out"; $new =~ s/OUT2/$out2/g;
	my $file1= $blast . "/" . $id; $new =~ s/CPU/$cpu/; $new =~ s/DB1/$file1/g; $new =~ s/DB/$id/g; $new =~ s/RAND/$random/g;  my $jobname = "$id-$random-RUN";
	$new =~ s/BLAST/tblastx/; $new =~ s/DIR/$out_dir/g;  $new =~ s/EMAIL/$email/; $new =~ s/JOBNAME/$jobname/;
	if (!$num)
	{
		my $out_id = $out_dir ."/" . $random . "_v_" . $id . ".tblastx.out";
		$new =~ s/OUT1/$out_id/;
		$new =~ s/TRINITY/$trinity/g;
		$new =~ s/NUM\.//g;
	}
	else
	{
		my $out_id = $out_dir ."/" . $random . "_v_" . $id . ".tblastx.\$SLURM_ARRAY_TASK_ID.out";
		$new =~ s/OUT1/$out_id/;
		my $tmp_trinity = $trinity  . ".\$SLURM_ARRAY_TASK_ID";
		$new =~ s/TRINITY/$tmp_trinity/g;
		$new =~ s/NUM/\$SLURM_ARRAY_TASK_ID/g; 
	}
	if (!-e $out_dir . "/" . $random . "_" . $id . ".start" && !-e $out_dir . "/" . $random . "_" . $id . ".finished")
	{  
		open(my $fo, ">", $out_dir . "/" .$random . "_v_" . $id. ".sh");
		print {$fo} $new;
		close($fo);
		my $tmp_sh = $out_dir . "/" .$random . "_v_" . $id . ".sh";
		`chmod a+x $tmp_sh`;
		if (!$test && !-e "$out_dir/". $random ."_".$id . ".sub") 
		{ 
			my $touch = "touch ". $random ."_".$id . ".sub";
			`$touch`;
			if (!$num) {`$cmd_run $tmp_sh`; }
			else {	
				`$cmd_multi $tmp_sh`; }
		}
	}
}

while ($ls_pblast =~ /([^\n\r]+)/g)
{
	my $file = $1; $file =~ s/\.phr//; $file =~ /([^\/]+)\Z/; my $id = $1; my $cpu2 = $cpu;
	if ($id eq "Avent") { $cpu2 = 10; }
	my $new = $in; my $file1= $blast . "/" . $id; print "$file1\n"; $new =~ s/DB1/$file1/g; $new =~ s/DB/$id/g;  $new =~ s/RAND/$random/g;
	$new =~ s/BLAST/blastx/; $new =~ s/DIR/$out_dir/; $new =~ s/EMAIL/$email/;  $new =~ s/CPU/$cpu2/;
my $out2 = $out_dir . "/" . $random . "_" . $id . ".out"; $new =~ s/OUT2/$out2/g;
	my $jobname = "$id-$random-RUN";
  	$new =~ s/JOBNAME/$jobname/;
	if (!$num)
	{
		my $out_id = $out_dir ."/" . $random . "_v_" . $id . ".blastx.out";
		$new =~ s/OUT1/$out_id/;
		$new =~ s/TRINITY/$trinity/g;
		$new =~ s/NUM\.//g;
	}
	else
	{
		my $out_id = $out_dir ."/" . $random . "_v_" . $id . ".tblastx.\$SLURM_ARRAY_TASK_ID.out";
		$new =~ s/OUT1/$out_id/;
		my $tmp_trinity = $trinity  . ".\$SLURM_ARRAY_TASK_ID";
		$new =~ s/TRINITY/$tmp_trinity/g;
 		$new =~ s/NUM/\$SLURM_ARRAY_TASK_ID/g; 
	}
        if (!-e $out_dir . "/" . $random . "_" . $id . ".start" && !-e $out_dir . "/" . $random . "_" . $id. ".finished")
	{
		open(my $fo, ">", $out_dir . "/" .$random . "_v_" . $id . ".sh");
		print {$fo} $new;
		close($fo);
		my $tmp_sh = $out_dir . "/" .$random . "_v_" . $id . ".sh";
		`chmod a+x $tmp_sh`;
		if (!$test && !-e "$out_dir/". $random ."_".$id . ".sub" && !-e "$out_dir/". $random ."_".$id . ".start" && !-e "$out_dir/". $random ."_".$id . ".finished")
		{
			my $touch = "touch $out_dir/". $random ."_".$id . ".sub";
                        `$touch`;
			if (!$num) { `$cmd_run $tmp_sh`; }
			else { `$cmd_multi $tmp_sh`; }
		}
	}
}

if (!$skip_nr)
{
my $new = $in3;
my $id = "NR";
my $jobname = "$id-$random-RUN";
$new =~ s/JOBNAME/$jobname/;
$new =~ s/RAND/$random/g;
$new =~ s/DIR/$out_dir/g;  $new =~ s/EMAIL/$email/; $new =~ s/CPU/$cpu/;
my $out2 = $out_dir . "/" . $random . "_" . $id . ".out"; $new =~ s/OUT2/$out2/g;
if (!$num)
{
        $new =~ s/TRINITY/$trinity/g;
 	$new =~ s/NUM//g;
}
else
{
	my $tmp_trinity = $trinity  . ".\$SLURM_ARRAY_TASK_ID";
        $new =~ s/TRINITY/$tmp_trinity/g;
	$new =~ s/NUM/\$SLURM_ARRAY_TASK_ID/g;
}
if (!-e $out_dir . "/" . $random . "_" . $id . ".start" && !-e $out_dir . "/" . $random . "_" . $id . ".finished")
{
	open(my $fo, ">", $out_dir . "/" .$random . "_v_" . $id. ".sh");
     print {$fo} $new;
     close($fo);
     my $tmp_sh = $out_dir . "/" .$random . "_v_" . $id . ".sh";
     `chmod a+x $tmp_sh`;
     if (!$test && !-e "$out_dir/". $random ."_".$id . ".sub" && !-e "$out_dir/". $random ."_".$id . ".start" && !-e "$out_dir/". $random ."_".$id . ".finished")
     {
     	my $touch = "touch $out_dir/". $random ."_" . $id. ".sub";
        `$touch`;
        if (!$num) {`$cmd_run $tmp_sh`; }
        else { `$cmd_multi $tmp_sh`; }
     }
}
}
if ($fq_dir)
{
	if (! -e $out_dir . "/" . $random .".rsemdb")
	{
		`module load rsem 
		 rsem-prepare-reference --bowtie $trinity $out_dir/$random`;
		`touch $out_dir/$random.rsemdb`;
	}  
  	if (! -e $out_dir . "/" . $random .".btdb")
        {
                `module load bowtie 
                 bowtie-build $trinity $out_dir/$random`;
                `touch $out_dir/$random.btdb`;
        }
	my $ls = `ls $fq_dir/*1.f*`;
	while ($ls =~ /([^\n\r+]+)/g)
	{
		my $fq = $1;
		$fq =~ s/$fq_dir//g;
		$fq =~ s/\A\///g;
		print $fq, "\n";
		if ($fq =~ /(.+)1\.f(.*)q(.*)/)
		{
			my $id = $1; my $end = "f" . $2 . "q" . $3;
			print  "$fq_dir/$id" . "2." . $end, "\n";
			if (-e "$fq_dir/$id" . "2." . $end)
			{
				my $new = $in2;
				my $jobname = "$id-$random-RUN";
				$new =~ s/JOBNAME/$jobname/;
				$new =~ s/RAND/$random/g; $new =~ s/CPU/$cpu/;
				$new =~ s/DIR/$out_dir/g;  $new =~ s/EMAIL/$email/;
				my $out1 = "$out_dir/$random" . "_$id" . "_rsem";
				my $id1 = "$fq_dir/$fq";
				my $id2 = "$fq_dir/$id" . "2.$end";
				my $db = "$out_dir/$random";
				my $out2 = $out_dir . "/" . $random . "_" . $id . ".out"; $new =~ s/OUT2/$out2/g;
				$new =~ s/OUT1/$out1/g; $new =~ s/ID1/$id1/g; $new =~ s/ID2/$id2/g; $new =~ s/DB1/$db/g; $new =~ s/ID/$id/g;
				if (!-e $out_dir . "/" . $random . "_" . $id . ".rsemstart" && !-e $out_dir . "/" . $random . "_" . $id . ".rsemdone")
				{
        				open(my $fo, ">", $out_dir . "/" .$random . "_v_" . $id. ".sh");
     					print {$fo} $new;
     					close($fo);
     					my $tmp_sh = $out_dir . "/" .$random . "_v_" . $id . ".sh";
     					`chmod a+x $tmp_sh`;
     					if (!$test && !-e "$out_dir/". $random ."_".$id . ".sub" && !-e "$out_dir/". $random ."_".$id . ".rsemstart" && !-e "$out_dir/". $random ."_".$id . ".rsemdone")
     					{
        					my $touch = "touch $out_dir/". $random ."_" . $id. ".sub";
        					`$touch`;
        					`$cmd_run $tmp_sh`;
     					}
				}

			}
		}
	}
}


