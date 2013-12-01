#!/usr/bin/perl

use strict;
use warnings;
use autodie;

my $N = shift || 20;
my $M = $N * $N;

while( <> ) {
    last if /satisfiable/i;
}
if( /unsatisfiable/i ) {
    print "NA\n";
    exit;
}

my( @F );
my $s = join " ", <>;
$s =~ s/(-?\d+)/process($1)/egs;

print join ",", @F[1..$M];
print "\n";

sub process {
    my( $x ) = @_;
    set( decode($x), $x > 0 );
}

sub set {
    my( $i, $k ) = @_;
    $F[$i] = 0+$k;
}

sub decode {
    my( $x ) = @_;
    return abs($x);
}
