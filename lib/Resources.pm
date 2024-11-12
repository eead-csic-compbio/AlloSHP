package Resources;

use strict;
require Exporter;

our @ISA = qw( Exporter );
our @EXPORT = qw( 
  get_time_RAM
);

# takes Memory::Usage report and return two strings: i) time in seconds, ii) RAM in Mb
sub get_time_RAM {
  my ($report) = @_;
  my ($time, $RAM);

  while($report =~ /(\d+)\s+(\d+)\s+\(\s*\d+\)\s+(\d+)\s+\(\s*\d+\)\s+(\d+)\s+\(\s*\d+\)\s+(\d+)\s+\(\s*\d+\)\s+(\d+)\s+\(\s*\d+\)/g) {
    $time = $1;
    $RAM = sprintf("%2.1f",($2+$3+$4+$5+$6)/1024);
  }

  return ($time, $RAM);
}
