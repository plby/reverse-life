#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use File::Copy;
use Time::HiRes qw(gettimeofday);

my( $directory, $jump, $MAX );
$directory = shift;
$jump = shift || 2;
$MAX = 50_000;

chdir $directory;

OUTER:
for my $i ( 1 .. $MAX ) {
	print "i=$i\n";
	chdir $i;
	next OUTER unless -e "delta" and -e "stop";
	my $delta = `cat delta`;
	print "delta=$delta\n";
	for my $j ( 0 .. $delta-1 ) {
		print "j=$j\n";
		my $prev = "frog$jump-$j";
		$prev = "stop" if $j == 0;
		my $next = "frog$jump-".($j+1);

		my $command = "cat $prev | perl ../../../src/undo_k.pl 20 $jump 1 | picosat | perl ../../../src/read_sat.pl > $next";
		print "$command\n";
		my($starts, $startus) = gettimeofday;
		system( $command );
		my($stops , $stopus ) = gettimeofday;
		my $time = ($stops-$starts) + ($stopus-$startus)/1_000_000;
		printf "Computation took %.2f seconds.\n", $time;

		my $data = `cat $next`;
		print $data;
		if( $data =~ m/NA/i ) {
			print "j=$j\tFAIL\n";
			next OUTER;
		} elsif( $data =~ m/[01]/ ) {
			print "j=$j\t\tSUCCESS\n";
		} else {
			warn "This shouldn't happen.  Debug me?\n";
			next OUTER;
		}
	}

	print "COMPLETE INVERSION up to delta=$delta!\n";
	copy( "frog$jump-$delta", "start-grid" ) unless -e "start-grid";
} continue {
	chdir "..";
}
