#!/usr/bin/perl -w
use strict;

my $TMP  = "/tmp/cgiknot-".getpwuid($<);
my $file = shift;

# my ($REC,$NUM,$ATM,$RES,$CHN,$NRES,$X,$Y,$Z) = (0,1,3,5,6,8,10,11,12);
# my ($MN1, $MN2, $ELM) = (13,14,2);

print "Content-Type: text/html; charset=ISO-8859-1\n\n";

open(my $fh, "<", "$TMP/$file") || die "Can not open file: $!";
my @atom = ();
while(<$fh>){
    my @t = split /\s+/;
    push @atom, \@t if @t == 4;
}
close $fh;

if(@atom>2){
    for(my $i=0; $i<@atom; $i++){
	my $t = $atom[$i];
	if($i == 0  || $i == $#atom){
	    my $j = $i==0 ? 1 : $i-1;
	    for(my $k=1;$k<4; $k++){
		$t->[$k] = $t->[$k] + ($atom[$j][$k]-$t->[$k])*0.6;
	    }
	}
    }

    my @tatom = ();
    for(my $i=0; $i<@atom-1; $i++){
	my ($d,$cosx,$cosy,$cosz) = atom_dist($atom[$i], $atom[$i+1]);

	my $t = $atom[$i];
	push @tatom,$t;
	if($d>3.8){
	    for(my $dj=3.768; $dj<$d; $dj+=3.768){
		push @tatom, [0, $t->[1]+$dj*$cosx, $t->[2]+$dj*$cosy, $t->[3]+$dj*$cosz];

	    }
	}
    }
    push @tatom,$atom[$#atom];
    @atom = @tatom;
}

for(my $i=0; $i<@atom; $i++){
    my $t = $atom[$i];

    printf("%-6s %4d %1s%-3s %3s %1s %3d     %7.3f %7.3f %7.3f  %4.2f%6.2f %11s\n", 'ATOM', $i, '', 'CA', 'GLY', 'A', $i, $t->[1]/50, $t->[2]/50, $t->[3]/50, 1.00, 99.00, 'C');
}

sub atom_dist{
    my ($a1, $a2)=@_;
    my ($dx,$dy,$dz);
    my $ret = sqrt(
	($dx=($a2->[1] - $a1->[1]))**2 +
	($dy=($a2->[2] - $a1->[2]))**2 +
	($dz=($a2->[3] - $a1->[3]))**2 
    );
    return ($ret, $dx/$ret,  $dy/$ret, $dz/$ret) if wantarray;
    return $ret;
}

