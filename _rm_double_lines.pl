#!/usr/bin/perl -w
use strict;

open(FH,"zcat $ARGV[0] |"); #){ die "# usage: $0 <file.vcf.gz>\n" }
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
