#!/usr/bin/perl

package ReverseLife;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = "1.00";
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(display_grid);
%EXPORT_TAGS = ( DEFAULT => [qw(&display_grid)],
	       );

sub display_grid {
	my %defaults = ( N => 20 );
	my( %args ) = @_;
	%args = (%defaults, %args );

	my $r = "";
	for my $i ( 0 .. $args{N}-1 ) {
		for my $j ( 0 .. $args{N}-1 ) {
			$r .= display_char( %args, c => $args{flat}->[$args{N}*$i + $j]);
		}
		$r .= "\n";
	}
	return $r;
}

sub display_char {
	my %defaults = ( N => 20, chars => "01" );
	my( %args ) = @_;
	%args = (%defaults, %args );

	return substr $args{chars}, $args{c}, 1;
}

1;
