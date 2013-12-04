#!/usr/bin/perl

use strict;
use warnings;
use autodie;

my( $directory, $MAX );
$directory = shift;
$MAX = 50_000;

chdir $directory;

my( $total, @tau );
@tau = (0) x 6;
OUTER:
for my $i ( 1 .. $MAX ) {
	chdir $i;
	next OUTER unless -e "delta" and -e "stop";
	my $delta = 0+`cat delta`;
	next OUTER unless $delta == 5;

	last OUTER unless -e "frog2-1";
	$total++;
	if( -e "start-grid" ) {
		$tau[5]++;
	} elsif( -e "frog2-5" ) {
		$tau[4]++;
	} elsif( -e "frog2-4" ) {
		$tau[3]++;
	} elsif( -e "frog2-3" ) {
		$tau[2]++;
	} elsif( -e "frog2-2" ) {
		$tau[1]++;
	} else {
		$tau[0]++;
	}
	last if $total == 100;
} continue {
	chdir "..";
}

print "$total\n";
for my $i ( 0 .. 5 ) {
	print "$i\t$tau[$i]\n";
}
