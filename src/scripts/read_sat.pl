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
    last if /satisfiable/i;
}
if( /unsatisfiable/i ) {
    print "NA\n";
    exit;
}

my( @F );
my $s = join " ", <>;
$s =~ s/(-?\d+)/process($1)/egs;

shift @F;
print display_grid(flat => \@F);

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
