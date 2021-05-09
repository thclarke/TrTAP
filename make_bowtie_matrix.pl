#!/usr/bin/perl
use Getopt::Long;
use strict;
use warnings;
use File::Basename;

my ($dir, $out, $add, $help, $index, $end, $paired, $col, $use_new);
    GetOptions("n|new"=>\$use_new, "d|dir:s" =>\$dir, "o|out:s"=>\$out, "a|add:s"=>\$add, "c|col:s"=>\$col, "i|index:s"=>\$index, "paired|p"=>\$paired, "e|end:s"=>\$end, "h|?|help"=>\$help);

if ($help || !($dir))
{
    print STDERR "Make RSEM Matrix\n";
    print STDERR "-----------------USAGE--------------\n";
    print STDERR " -d directory with bowtie results \n";
    print STDERR " -o output file\n";
    print STDERR " -c column in bowtie count file to use\n";
    print STDERR " -n use new counts\n";
    print STDERR " -a addition to the start of the file. Removed for output column name\n";
    die();
}

my $config_file = dirname($ARGV[0]) . "/trtap.ini";
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
my $new;
my $a = $add;
my $list;
opendir(DIR, $dir);
my @f = readdir(DIR);
my $lst;
foreach my $b (@f)
{
    my $tmp;
    if ($a)
    {
    	if ($b =~ /$a(\S+)$end.bowtie([s2]*).cnt(s*)\Z/)
    	{
       	 $tmp = $1;
       	 $lst->[scalar(@$lst)] = $tmp; 
        }
    }
    else
   {
	if ($use_new)
	{
		if ($b =~ /(\S+)$end.new_bowie.cnt(s*)2\Z/)
        	{
         	   $tmp = $1;
         	   $lst->[scalar(@$lst)] = $tmp;
        	}	
	}
	else {
	if ($b =~ /(\S+)$end.bowtie([s2]*).cnt(s*)\Z/)
        {
            $tmp = $1;
           $lst->[scalar(@$lst)] = $tmp;
        } }
   }
   if ($tmp)
   {
	print STDERR "$b\n";
	open(FILE, "<", $dir ."/" .$b) or die();
        <FILE>;
        while(!eof(FILE))
        {
                my $v = <FILE>; chomp($v);
                my @n = split "\t", $v;
                #if ($n[0] =~ /\A(\d)/) { if ($add) { $n[0] = $add ."_". $n[0]; }}
                $list->{$n[0]}->{$tmp} = $n[$col];
        }
    }
}
open(my $fo, ">", $out);
print {$fo} "Name";
for (my $i = 0; $i < scalar(@$lst); $i++)
{
        print {$fo} "\t", $lst->[$i];
}
print {$fo} "\n";
foreach my $a (keys(%$list))
{
    
                print {$fo} $a;
                for (my $i = 0; $i < scalar(@$lst); $i++)
                {
			if ($list->{$a}->{$lst->[$i]})
			{
                        	print {$fo} "\t", 0+$list->{$a}->{$lst->[$i]};
			}
			else 
			{
				print {$fo} "\t0";
			}
                }
                print {$fo} "\n";
     
}
