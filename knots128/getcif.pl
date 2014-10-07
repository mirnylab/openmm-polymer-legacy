#!/usr/bin/perl -w
use strict;

my $TMP  = "/tmp/cgiknot-".getpwuid($<);
my $file = shift;

my ($REC,$NUM,$ATM,$RES,$CHN,$NRES,$X,$Y,$Z) = (0,1,3,5,6,8,10,11,12);
my ($MN1, $MN2, $ELM) = (13,14,2);

print "Content-Type: text/html; charset=ISO-8859-1\n\n";

my $pdb = 1;
# $pdb = 1 if $file=~/\.pdb/;

open(my $fh, "<", "$TMP/$file") || die "Can not open file: $!";
my $atom=0;
while(<$fh>){
    if($pdb){
	print; next;
    }
    if(/^\s*_atom/){ $atom++; next };
    if($atom && /^(ATOM|HETATM)/){
	my @t = split /\s+/;
	my $xatm = " ";
	$t[$ATM] =~ s/\'//g;
	if(length $t[$ATM]>3){
	    $xatm = substr($t[$ATM],3,1);
	    $t[$ATM] = substr($t[$ATM],0,3);
	}

	printf("%-6s%5d %1s%-3s %3s %1s%4s     %7.3f %7.3f %7.3f  %4.2f%6.2f %11s\n",
	    $t[$REC], $t[$NUM], $xatm, $t[$ATM], $t[$RES], $t[$CHN], $t[$NRES],
	    $t[$X], $t[$Y], $t[$Z], $t[$MN1], $t[$MN2], $t[$ELM]
	);
    }
}
close $fh;
