#!/usr/bin/perl

# Read test.csv or train.csv and generate one directory for each
# position in them

use strict;
use warnings;
use autodie;

use File::Basename;
use lib dirname (__FILE__);
use ReverseLife qw(display_grid);

my $N = shift || 20;
my $M = $N * $N;

my( $header, $trainq );
chomp( $header = <> );
$trainq = $header =~ /stop/;

while( <> ) {
	chomp;
	my( @F, $id, $delta, @start, @stop );
	@F = split ",", $_;
	$id    = shift @F;
	$delta = shift @F;
	@stop  = splice @F, -$M if $trainq;
	@start = @F;

	mkdir $id unless -e $id;
	open DELTA, ">", "$id/delta";
	print DELTA $delta, "\n";
	close DELTA;

	open START, ">", "$id/start";
	print START display_grid(flat => \@start);
	close START;

	if( $trainq ) {
		open STOP, ">", "$id/stop";
		print STOP display_grid(flat => \@stop);
		close STOP;
	}
}
