#!/usr/bin/perl

# This file accepts a CNF-SAT input and then determines whether any of
# the variable values are forced.  It assumes the input begins with
# the "p cnf VARS CLAUSES" line.

use File::Temp qw(tempfile);

use strict;
use warnings;
use autodie;

# Process arguments
my( $solver, @args, $min_var, $max_var );
$solver = shift || "picosat";
$min_var = shift;
$max_var = shift;

push @args, "-model" if $solver =~ m/minisat|glucose/;

my( $header, $VARS, $CLAUSES, $REST );
# Read header
chomp( $header = <> );
($VARS, $CLAUSES) = $header =~ m/p cnf (\d+) (\d+)/;
$min_var ||= 1;
$max_var ||= $VARS;

# Read the clauses
$REST = join "", <>;

# Initialize a couple temporary files for input and output
my( $infile, $outfile );
(undef, $infile ) = tempfile( "possXXXX", DIR => "/tmp", UNLINK => 0, SUFFIX => ".sat.in"  );
(undef, $outfile) = tempfile( "possXXXX", DIR => "/tmp", UNLINK => 0, SUFFIX => ".sat.out" );
#warn "$infile $outfile\n";

# Write back SAT problem
{
	open IN, ">", $infile;
	print IN "p cnf $VARS $CLAUSES\n";
	print IN $REST;
	close IN;
}

# Find one possible value of the variables
solve();

# Read assignment
my( @one );
my $sat = read_assignment( $outfile, \@one );
die "Expected solvable SAT instance." unless $sat;

# Loop to find remaining possibilities
my( @opposite_possible );
for my $i ( $min_var .. $max_var ) {
	$opposite_possible[$i] = 0;
}

while( 1 ) {
	open IN, ">", $infile;
	print IN "p cnf $VARS ", (1+$CLAUSES), "\n";
	for my $i ( $min_var .. $max_var ) {
		next if $opposite_possible[$i];
		print IN "-" if $one[$i];
		print IN "$i ";
	}
	print IN "0\n";
	print IN $REST;
	close IN;

	solve();

	my( @two );
	my $sat = read_assignment( $outfile, \@two );
	last unless $sat;

	for my $i ( $min_var .. $max_var ) {
		next unless $two[$i] != $one[$i];
		$opposite_possible[$i] = 1;
	}
}

unlink $infile;
unlink $outfile;

# Output forced variables
for my $i ( $min_var .. $max_var ) {
	next if $opposite_possible[$i];
	print "-" unless $one[$i];
	print "$i\n";
}

sub read_assignment {
	my( $file, $array ) = @_;
	my( $result );

	open OUT, "<", $file;
	my $response = join " ", <OUT>;
	close OUT;
	if( $response =~ m/unsatisfiable/i ) {
		$result = 0;
	} elsif( $response =~ m/satisfiable/i ) {
		$result = 1;
	} else {
		print $response;
		die "Bad SAT response.\n";
	}
	$response =~ s/(-?\d+)/$array->[abs($1)] = ($1 > 0 ? 1 : 0)/ge;
	return $result;
}

my( $count );
sub solve {
	$count++;
	warn "Solving instance $count.\n";
	system( $solver, @args, "-o", $outfile, $infile );
}
