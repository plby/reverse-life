#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use File::Basename;
use lib dirname (__FILE__);
use ReverseLife qw(display_grid);

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
			print display_grid( N     => $N,
					    flat  => [splice @F, 0, $M],
					    chars => ".#"
					  );
			print "\n";
		}
	}
	print "---\n";
}
