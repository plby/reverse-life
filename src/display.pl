#!/usr/bin/perl

use strict;
use warnings;
use autodie;

my $N = shift || 20;
my $M = $N * $N;

while( <> ) {
	chomp;

	if( /NA/i ) {
		print "NA\n";
	} else {
		my( @F );
		@F = split ",", $_;
		splice @F, 0, @F - $M * int(@F/$M);
		while( @F ) {
			print display( splice @F, 0, $M );
		}
	}
	print "---\n";
}

sub display {
	my( @x ) = @_;

	my $r = "";
	for my $i ( 0 .. $N-1 ) {
		for my $j ( 0 .. $N-1 ) {
			$r .= char($x[$N*$i + $j]);
		}
		$r .= "\n";
	}
	$r .= "\n";
	return $r;
}

sub char {
	return substr ".#", shift, 1;
}
