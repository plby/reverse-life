#!/usr/bin/perl

# This file accepts a CNF-SAT input and then determines whether any of
# the variable values are forced.  It assumes the input begins with
# the "p cnf VARS CLAUSES" line.

use File::Temp qw(tempfile);

use strict;
use warnings;
use autodie;

my( $solver, $header, $VARS, $CLAUSES, $REST );

$solver = shift || "picosat";

# Read header
chomp( $header = <> );
($VARS, $CLAUSES) = $header =~ m/p cnf (\d+) (\d+)/;

# Read the clauses
$REST = join "\n", <>;

# Initialize a couple temporary files for input and output
my( $infile, $outfile );
(undef, $infile ) = tempfile( "possXXXX", DIR => "/tmp", UNLINK => 0, SUFFIX => ".sat.in"  );
(undef, $outfile) = tempfile( "possXXXX", DIR => "/tmp", UNLINK => 0, SUFFIX => ".sat.out" );

# Find one possible value of the variables
system( $solver, "-o", $outfile, $infile );

# 
