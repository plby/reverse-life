#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use File::Copy;
use Time::HiRes qw(gettimeofday);

my( $directory, $limit );
$directory = shift;
$limit = shift || 5;

chdir $directory;

OUTER:
for my $i ( "00001" .. "50000" ) {
	my $subdir = join "/", split "", $i;
	chdir $subdir;

	next OUTER unless -e "delta" and -e "stop";
	next OUTER if -e "start-grid" and (-s "start-grid" >= 50);
	my $delta = 0+`cat delta`;
	next OUTER unless 1 <= $delta and $delta <= 5;
	next OUTER if $delta > $limit;

	my $prefix = "~/i/reverse-life/src";
	my $command = "cat stop | perl $prefix/undo_k.pl 20 $delta | picosat | perl $prefix/read_sat.pl > start-grid";
	my($starts, $startus) = gettimeofday;
	system( $command );
	my($stops , $stopus ) = gettimeofday;
	my $time = ($stops-$starts) + ($stopus-$startus)/1_000_000;
	printf "$i\t$delta\t%.2f\n", $time;
} continue {
	chdir "../../../../..";
}
