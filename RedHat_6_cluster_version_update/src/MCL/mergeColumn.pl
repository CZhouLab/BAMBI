#! /usr/bin/perl -w

#Function: Merge two files into one
#Input: file1 and file2, where each file is formated as columns 
#       and 
#       two files have the same rows.
#Output: (to stdout by default) columns of file1 and file2, where columns are delimited by tab.
#
#Author:ZHOU Chan
#Email:zhouchan99@gmail.com
#Date:2007-06-03

use strict;
my $usage="$0 file1 file2 > outfile\n";

if(@ARGV<2){die $usage;}
open IF1, $ARGV[0];
open IF2, $ARGV[1];
chomp(my @f1lines=<IF1>);
chomp(my @f2lines=<IF2>);
my $f1linenumber=@f1lines;
my $f2linenumber=@f2lines;
if($f1linenumber ne $f2linenumber){die "Requirement: file1 and file2 should have the same lines\n";}

for(my $i=0;$i<$f1linenumber;$i++){
      print "$f1lines[$i]\t$f2lines[$i]\n";
}

close IF1;
close IF2;

