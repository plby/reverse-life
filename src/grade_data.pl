#!/usr/bin/perl

use strict;
use warnings;
use autodie;

my( $train, $submission, $N ) = @ARGV;
die "Usage: $0 TRAIN_FILE SUBMISSION_FILE [N]\n"
  unless defined $train and defined $submission
  and -e $train and -e $submission;
$N ||= 20;
my $M = $N*$N;

open TRAIN     , "<", $train     ;
open SUBMISSION, "<", $submission;

# Check header
{
	<TRAIN>;
	my $t = <DATA>;
	my $s = <SUBMISSION>;
	die "Bad header in submission file.\n" unless $s eq $t;
}

my( $right, $total );
for my $i ( 1 .. 50_000 ) {
	my( @s, @t );
	{
		my $t = <TRAIN>;
		my $s = <SUBMISSION>;
		@t = split ",", $t;
		@s = split ",", $s;
	}
	splice @t, 1, 1; # the training data includes a delta in it
	die "ID numbers don't match up.\n" unless (splice @t, 0, 1) == (splice @s, 0, 1);

	# Grade!
	for my $j ( 1 .. $M ) {
		my $t = shift @t;
		my $s = shift @s;
		$right++ if $s == $t;
		$total++;
	}
}
die "Didn't find any cells to grade!?\n"
  unless $total > 0;

my $result = $right / $total;
print "$result\n";

__DATA__
id,delta,start.1,start.2,start.3,start.4,start.5,start.6,start.7,start.8,start.9,start.10,start.11,start.12,start.13,start.14,start.15,start.16,start.17,start.18,start.19,start.20,start.21,start.22,start.23,start.24,start.25,start.26,start.27,start.28,start.29,start.30,start.31,start.32,start.33,start.34,start.35,start.36,start.37,start.38,start.39,start.40,start.41,start.42,start.43,start.44,start.45,start.46,start.47,start.48,start.49,start.50,start.51,start.52,start.53,start.54,start.55,start.56,start.57,start.58,start.59,start.60,start.61,start.62,start.63,start.64,start.65,start.66,start.67,start.68,start.69,start.70,start.71,start.72,start.73,start.74,start.75,start.76,start.77,start.78,start.79,start.80,start.81,start.82,start.83,start.84,start.85,start.86,start.87,start.88,start.89,start.90,start.91,start.92,start.93,start.94,start.95,start.96,start.97,start.98,start.99,start.100,start.101,start.102,start.103,start.104,start.105,start.106,start.107,start.108,start.109,start.110,start.111,start.112,start.113,start.114,start.115,start.116,start.117,start.118,start.119,start.120,start.121,start.122,start.123,start.124,start.125,start.126,start.127,start.128,start.129,start.130,start.131,start.132,start.133,start.134,start.135,start.136,start.137,start.138,start.139,start.140,start.141,start.142,start.143,start.144,start.145,start.146,start.147,start.148,start.149,start.150,start.151,start.152,start.153,start.154,start.155,start.156,start.157,start.158,start.159,start.160,start.161,start.162,start.163,start.164,start.165,start.166,start.167,start.168,start.169,start.170,start.171,start.172,start.173,start.174,start.175,start.176,start.177,start.178,start.179,start.180,start.181,start.182,start.183,start.184,start.185,start.186,start.187,start.188,start.189,start.190,start.191,start.192,start.193,start.194,start.195,start.196,start.197,start.198,start.199,start.200,start.201,start.202,start.203,start.204,start.205,start.206,start.207,start.208,start.209,start.210,start.211,start.212,start.213,start.214,start.215,start.216,start.217,start.218,start.219,start.220,start.221,start.222,start.223,start.224,start.225,start.226,start.227,start.228,start.229,start.230,start.231,start.232,start.233,start.234,start.235,start.236,start.237,start.238,start.239,start.240,start.241,start.242,start.243,start.244,start.245,start.246,start.247,start.248,start.249,start.250,start.251,start.252,start.253,start.254,start.255,start.256,start.257,start.258,start.259,start.260,start.261,start.262,start.263,start.264,start.265,start.266,start.267,start.268,start.269,start.270,start.271,start.272,start.273,start.274,start.275,start.276,start.277,start.278,start.279,start.280,start.281,start.282,start.283,start.284,start.285,start.286,start.287,start.288,start.289,start.290,start.291,start.292,start.293,start.294,start.295,start.296,start.297,start.298,start.299,start.300,start.301,start.302,start.303,start.304,start.305,start.306,start.307,start.308,start.309,start.310,start.311,start.312,start.313,start.314,start.315,start.316,start.317,start.318,start.319,start.320,start.321,start.322,start.323,start.324,start.325,start.326,start.327,start.328,start.329,start.330,start.331,start.332,start.333,start.334,start.335,start.336,start.337,start.338,start.339,start.340,start.341,start.342,start.343,start.344,start.345,start.346,start.347,start.348,start.349,start.350,start.351,start.352,start.353,start.354,start.355,start.356,start.357,start.358,start.359,start.360,start.361,start.362,start.363,start.364,start.365,start.366,start.367,start.368,start.369,start.370,start.371,start.372,start.373,start.374,start.375,start.376,start.377,start.378,start.379,start.380,start.381,start.382,start.383,start.384,start.385,start.386,start.387,start.388,start.389,start.390,start.391,start.392,start.393,start.394,start.395,start.396,start.397,start.398,start.399,start.400,stop.1,stop.2,stop.3,stop.4,stop.5,stop.6,stop.7,stop.8,stop.9,stop.10,stop.11,stop.12,stop.13,stop.14,stop.15,stop.16,stop.17,stop.18,stop.19,stop.20,stop.21,stop.22,stop.23,stop.24,stop.25,stop.26,stop.27,stop.28,stop.29,stop.30,stop.31,stop.32,stop.33,stop.34,stop.35,stop.36,stop.37,stop.38,stop.39,stop.40,stop.41,stop.42,stop.43,stop.44,stop.45,stop.46,stop.47,stop.48,stop.49,stop.50,stop.51,stop.52,stop.53,stop.54,stop.55,stop.56,stop.57,stop.58,stop.59,stop.60,stop.61,stop.62,stop.63,stop.64,stop.65,stop.66,stop.67,stop.68,stop.69,stop.70,stop.71,stop.72,stop.73,stop.74,stop.75,stop.76,stop.77,stop.78,stop.79,stop.80,stop.81,stop.82,stop.83,stop.84,stop.85,stop.86,stop.87,stop.88,stop.89,stop.90,stop.91,stop.92,stop.93,stop.94,stop.95,stop.96,stop.97,stop.98,stop.99,stop.100,stop.101,stop.102,stop.103,stop.104,stop.105,stop.106,stop.107,stop.108,stop.109,stop.110,stop.111,stop.112,stop.113,stop.114,stop.115,stop.116,stop.117,stop.118,stop.119,stop.120,stop.121,stop.122,stop.123,stop.124,stop.125,stop.126,stop.127,stop.128,stop.129,stop.130,stop.131,stop.132,stop.133,stop.134,stop.135,stop.136,stop.137,stop.138,stop.139,stop.140,stop.141,stop.142,stop.143,stop.144,stop.145,stop.146,stop.147,stop.148,stop.149,stop.150,stop.151,stop.152,stop.153,stop.154,stop.155,stop.156,stop.157,stop.158,stop.159,stop.160,stop.161,stop.162,stop.163,stop.164,stop.165,stop.166,stop.167,stop.168,stop.169,stop.170,stop.171,stop.172,stop.173,stop.174,stop.175,stop.176,stop.177,stop.178,stop.179,stop.180,stop.181,stop.182,stop.183,stop.184,stop.185,stop.186,stop.187,stop.188,stop.189,stop.190,stop.191,stop.192,stop.193,stop.194,stop.195,stop.196,stop.197,stop.198,stop.199,stop.200,stop.201,stop.202,stop.203,stop.204,stop.205,stop.206,stop.207,stop.208,stop.209,stop.210,stop.211,stop.212,stop.213,stop.214,stop.215,stop.216,stop.217,stop.218,stop.219,stop.220,stop.221,stop.222,stop.223,stop.224,stop.225,stop.226,stop.227,stop.228,stop.229,stop.230,stop.231,stop.232,stop.233,stop.234,stop.235,stop.236,stop.237,stop.238,stop.239,stop.240,stop.241,stop.242,stop.243,stop.244,stop.245,stop.246,stop.247,stop.248,stop.249,stop.250,stop.251,stop.252,stop.253,stop.254,stop.255,stop.256,stop.257,stop.258,stop.259,stop.260,stop.261,stop.262,stop.263,stop.264,stop.265,stop.266,stop.267,stop.268,stop.269,stop.270,stop.271,stop.272,stop.273,stop.274,stop.275,stop.276,stop.277,stop.278,stop.279,stop.280,stop.281,stop.282,stop.283,stop.284,stop.285,stop.286,stop.287,stop.288,stop.289,stop.290,stop.291,stop.292,stop.293,stop.294,stop.295,stop.296,stop.297,stop.298,stop.299,stop.300,stop.301,stop.302,stop.303,stop.304,stop.305,stop.306,stop.307,stop.308,stop.309,stop.310,stop.311,stop.312,stop.313,stop.314,stop.315,stop.316,stop.317,stop.318,stop.319,stop.320,stop.321,stop.322,stop.323,stop.324,stop.325,stop.326,stop.327,stop.328,stop.329,stop.330,stop.331,stop.332,stop.333,stop.334,stop.335,stop.336,stop.337,stop.338,stop.339,stop.340,stop.341,stop.342,stop.343,stop.344,stop.345,stop.346,stop.347,stop.348,stop.349,stop.350,stop.351,stop.352,stop.353,stop.354,stop.355,stop.356,stop.357,stop.358,stop.359,stop.360,stop.361,stop.362,stop.363,stop.364,stop.365,stop.366,stop.367,stop.368,stop.369,stop.370,stop.371,stop.372,stop.373,stop.374,stop.375,stop.376,stop.377,stop.378,stop.379,stop.380,stop.381,stop.382,stop.383,stop.384,stop.385,stop.386,stop.387,stop.388,stop.389,stop.390,stop.391,stop.392,stop.393,stop.394,stop.395,stop.396,stop.397,stop.398,stop.399,stop.400
