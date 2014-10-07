#!/usr/bin/perl -w
use strict;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use lib $ENV{KGB_MODULES};
use KGBCIF;
use POSIX qw/ INT_MIN/;

my $TMP = "/tmp";

my $iscif = 0;
my $straight = 0;
my $redknotbase = undef;
my $xyz;

GetOptions(
    "c!" => \$iscif,
    "s!" => \$straight,
    "x!" => \$xyz,
    "r=s"=> \$redknotbase,
);

my $p; 

if($iscif){
    $p = get_cif(\*ARGV);
}
elsif($xyz){
    $p = get_xyz(\*ARGV);
}
else{
    $p = get_pdb(\*ARGV);
}


foreach my $chn (sort {$a cmp $b} keys %$p){

    if($straight){
	run_knot($p->{$chn},$chn);
    }
    else{
	my @t=($p->{$chn}[0]);
	for(my $i=1; $i<@{$p->{$chn}}; $i++){
	    if($p->{$chn}[$i]{n} - $p->{$chn}[$i-1]{n} > 1 &&
		atom_dist($p->{$chn}[$i], $p->{$chn}[$i-1])>4 ){
		run_knot(\@t,$chn) if @t>1;
		#@t=($p->{$chn}[$i]{n});
		@t = ();
		# next;
	    }

	    push @t, $p->{$chn}[$i];
	}
	run_knot(\@t,$chn) if @t>1;
    }

}

sub run_knot{
    my ($chn, $chnsym) = @_;
    my ($conf_fh, $conf_file) = tempfile(DIR=>$TMP, SUFFIX=> ".conf");


    my $pmap = {};
    {
	my $i=1;
	foreach my $ca (@$chn){
	    $pmap->{$i++} = $ca->{n};
	}
    }

    $chnsym = "" if $chnsym =~ /^\s*$/;



    print $conf_fh "t=0\n\n",scalar @$chn,"\n";

    my $i = 1;
    foreach my $ca (@$chn){
	print $conf_fh $i++," $ca->{x} $ca->{y} $ca->{z}\n";
    }
    close $conf_fh;

    print "'$chnsym' $chn->[0]{n} ".$chn->[$#$chn]{n}." ";

    my @arg = $conf_file;
    @arg = ("-r","$redknotbase.$chnsym-$chn->[0]{n}",@arg)
	if defined $redknotbase;
    my @r = split /\s+/,`./Reduce_knot20 @arg`;


    my ($tobeg,$toend) = ($r[2]-1,@$chn-$r[3]);

    unless($r[2] == 0 && $r[3] == 0){
	$r[2] = $pmap->{$r[2]};
	$r[3] = $pmap->{$r[3]};
    }


    my $minflanc = $tobeg<$toend ? $tobeg : $toend;
    print "@r $minflanc\n";


    unlink $conf_file;
}




sub get_pdb{
    my ($fh)  = @_;
    my $h={};
    my $model;

    my ($prev_chn, $prev_nres) = ('####',INT_MIN); 
    while(<$fh>){
	if(/^ATOM|^HETATM/){
	    if(/^(\S{1,6})\s*(\d+)\s+(\S{1,4})\s*(\S+) (.)\s*(\S+)\s*?(.{8})(.{8})(.{8})\s+/){
		my ($rec, $num, $atm, $res, $chn, $nres, $x, $y, $z) = 
		   ($1,   $2,   $3,   $4,   $5,   $6,    $7, $8, $9);

		next unless $atm eq 'CA';
		next if $chn eq $prev_chn && $nres == $prev_nres;

		map {s/^\s+//} ($x,$y,$z);

		$prev_chn  = $chn;
		$prev_nres = $nres;

		push @{$h->{$chn}}, {n=>$nres,x=>$x,y=>$y,z=>$z};
	    }
	    next;
	}

	/^ENDMDL/ && ++$model && last; # first model is analyzed for now.
    }

    return $h;
}

sub get_cif{
    my ($fh) = @_;

    my $cif = read_cif($fh,cacheoff => 1);
    my $h   = {};

    foreach my $p (values %{$cif->{chain}}){
	next unless defined $p->{type} && $p->{type} =~ /polypeptide/i;
	foreach my $chn (keys %{$p->{pdbchain_model}}){
	    my $m = (sort {$a<=>$b} keys %{$p->{pdbchain_model}{$chn}})[0];
	    my $pch = $p->{pdbchain_model}{$chn}{$m};

	    my @ca  = ();
	    while( my ($resn, $res) = each %$pch){
		my $res_pdb  = $p->{pdbmap}{num}{$resn};
		next unless defined $res_pdb;

		my $resn_pdb;
		$resn_pdb = $1 if $res_pdb =~ /(\-?\d+)/;
		foreach my $atom (@$res){
		    if($atom->{sym} eq 'C'){
			push @ca,
			    {n=>$resn_pdb, x=>$atom->{x}, y=>$atom->{y}, z=>$atom->{z}};
			last;
		    }
		}
	    }
	    $h->{$chn} = [sort {$a->{n} <=> $b->{n}} @ca];
	}
    }
    return $h;
}

sub get_xyz{
    my ($fh) = @_;

    my $i = 1;
    my $h = {};

    while(<$fh>){
	my @t = split;
	next unless @t == 3;
	my ($x , $y, $z ) = split;
	push @{$h->{A}}, {n=> $i++, x=>$x, y=>$y, z=>$z};
    }
    return $h;

}

sub atom_dist{
    my ($a1, $a2) = @_;
    return sqrt(($a1->{x} - $a2->{x})**2 +($a1->{y} - $a2->{y})**2 +($a1->{z} - $a2->{z})**2);
}
