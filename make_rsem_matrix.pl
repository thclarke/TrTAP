#!/usr/bin/perl
use Getopt::Long;
use strict;
use warnings;

my ($dir, $out, $add, $help, $index, $end, $paired, $col);
$end = "";
GetOptions("d|dir:s" =>\$dir, "e|end:s"=>\$end, "o|out:s"=>\$out, "a|add:s"=>\$add, "c|col:s"=>\$col, "i|index:s"=>\$index, "paired|p"=>\$paired, "h|?|help"=>\$help);

if ($help || !($dir))
{
    print STDERR "Make RSEM Matrix\n";
    print STDERR "-----------------USAGE--------------\n";
    print STDERR " -d directory with rsem results \n";
    print STDERR " -o output file\n";
    print STDERR " -c column in rsem file to use\n";
    print STDERR " -e end rsem file to use\n";
    print STDERR " -a addition to the start of the file. Removed for output column name\n";
    die();
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
    	if ($b =~ /$a(\S+)$end.genes.results\Z/)
    	{
       	 $tmp = $1;
       	 $lst->[scalar(@$lst)] = $tmp; 
        }
    }
    else
   {
	if ($b =~ /(\S+)$end.genes.results\Z/)
        {
            $tmp = $1;
           $lst->[scalar(@$lst)] = $tmp;
        }
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
my $fo = *STDERR;
if ($out){
	open($fo, ">", $out);
}
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
                        print {$fo} "\t", 0+$list->{$a}->{$lst->[$i]};
                }
                print {$fo} "\n";
     
}
