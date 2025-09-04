#!/usr/bin/env perl
use strict;
use warnings;

# removes repeated VCF lines with indels, if any

# Input1: VCF file, which might be GZIP/BZIP2-compressed

# Bruno Contreras, Ruben Sancho EEAD-CSIC 2017

my $infile;
if($ARGV[0]){ $infile = $ARGV[0] }
else{ die "# usage: $0 <file.vcf[.gz|.bz2]>\n" }

# check input file format and open it accordingly
my $magic;
open(INFILE,$infile) || die "# cannot read $infile, exit\n";
sysread(INFILE,$magic,2);
close(INFILE);

if($infile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
{
  if(!open(FH,"gzip -dc $infile |"))
  {
    die "# cannot read GZIP compressed $infile $!\n"
      ."# please check gzip is installed\n";
  }
}
elsif($infile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
{
  if(!open(FH,"bzip2 -dc $infile |"))
  {
    die "# cannot read BZIP2 compressed $infile $!\n"
      ."# please check bzip2 is installed\n";
  }
}
else{ open(FH,"<",$infile) || die "# cannot read $infile $!\n"; }


my $prev_line = '';
while(my $this_line = <FH>) 
{
   if(!$prev_line) # first line ever or lines after skipped duplicates
   {
      $prev_line = $this_line;
   }
   else # remaining lines
   {
      my @this_data = split(/\t/,$this_line,3);
      my @prev_data = split(/\t/,$prev_line,3);

      if($#this_data == $#this_data)
      {
         if($this_data[0] eq $prev_data[0] && $this_data[1] == $prev_data[1])
         {
            $prev_line = '';
         }
         else
         {
            print $prev_line;
            $prev_line = $this_line;
         } 	 	
      }
      else
      {
        print $prev_line;
        $prev_line = $this_line;
      }
   }
}
close (FH);

print $prev_line;
