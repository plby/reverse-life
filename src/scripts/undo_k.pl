#!/usr/bin/perl

# This file generates a SAT instance that encodes a backwards life
# evolution of $K steps on an $N x $N grid.  If $flip is set, then it
# flips the order of the boards (which helps with other programs).

# Usage: $0 $N $K [$flip]

# Accepts, on standard input, the grid.  Everything but 0s and 1s is
# ignored, so the input is quite flexible.

# Example: echo "0100 0100 0100 0000" | perl undo_k.pl 4 2

use strict;
use warnings;
use autodie;

my( $N, $K, $flip );
$N = shift || 20;
$K = shift || 1;
$flip = shift;

my $M = $N * $N;

# Read the final position into a two-dimensional array @W
my( @W );
{
	$_ = join "", <>;
	my( @F );
	@F = m/([01])/gs;
	splice @F, 0, -$M;
	die "Failed to parse input.\n" unless @F == $M;
	for my $i ( 0 .. $N-1 ) {
		for my $j ( 0 .. $N-1 ) {
			$W[$i][$j] = shift @F;
		}
	}
}

# $K grids, each with $M = $N*$N variables
my( @clauses );
my $VARS = $M * $K;

# For convenience, we include a variable encoding true whose negation
# of course encodes false
$VARS++;
my $alive = $VARS;
add_clause( $alive );

for my $k ( 0 .. $K-1 ) {
	for my $i ( 0 .. $N-1 ) {
		for my $j ( 0 .. $N-1 ) {
			my( $next, @prev );
			$next = get_var( $k+1, $i, $j );
			for my $di ( -1 .. 1 ) {
				for my $dj ( -1 .. 1 ) {
					push @prev, get_var( $k, $i+$di, $j+$dj );
				}
			}
			encode_sat_step( $next, @prev );
		}
	}
}

output($VARS, @clauses);

# Subroutines related to the manipulation of the CNF-SAT structure
sub add_clause {
	my( @x ) = @_;
	push @clauses, [@x];
}

sub output {
	my( $VARS, @clauses ) = @_;
	print "p cnf $VARS ", 0+@clauses, "\n";
	for my $c ( @clauses ) {
		print "@$c 0\n";
	}
}

# Get the variable underlying a position on one of the (K+1) grids,
# where one of the grids is the final configuration.  Returns 0 if the
# position is off the board.  (Recall that 0 is an invalid variable
# number.)
sub get_var {
	my( $k, $i, $j ) = @_;
	return 0 if $i < 0 or $i >= $N or $j < 0 or $j >= $N;
	return ($W[$i][$j] ? $alive : -$alive) if $k == $K;
	$k = $K-1 - $k if $flip;
	return $k*$M + $N*$i + $j + 1;
}

# Encode a single step of the Life evolution as a bunch of SAT
# clauses.  Takes one successor location and nine predecessors in
# row-major order, with the predecessor center in the middle.  A zero
# represents off the grid and is handled natively.
sub encode_sat_step {
	my( $next, @prev ) = @_;
      OUTER:
	for my $i ( 0 .. 2**9-1 ) {
		my( @boolprev, $boolnext, @clause );
		@boolprev = map {($i & (1 << $_)) ? 1 : 0}  reverse 0 .. 8;
		for my $j ( 0 .. 8 ) {
			next OUTER if $prev[$j] == 0 and $boolprev[$j] == 1;
		}
		$boolnext = bool_step(@boolprev);

		push @clause, ($boolnext ? 1 : -1) * $next;
		for my $j ( 0 .. 8 ) {
			push @clause, ($boolprev[$j] ? -1 : 1) * $prev[$j]
			  unless $prev[$j] == 0;
		}
		add_clause( @clause );
	}
}

# A single step of the Life evolution, as implementing with Boolean
# variables.  Same arguments as the previous routine, except
# everything is 0/1 here.
sub bool_step {
	my( @prev ) = @_;
	my $center = splice @prev, 4, 1;
	my $neighbors = 0;
	for my $p ( @prev ) {
		$neighbors += $p;
	}
	return $center if $neighbors == 2;
	return 1       if $neighbors == 3;
	return 0;
}
