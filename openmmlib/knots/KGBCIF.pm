# $Id: KGBCIF.pm,v 1.5 2007/02/15 20:53:46 kolesov Exp $
package KGBCIF;
require Exporter;
@ISA    = qw(Exporter);
@EXPORT = qw(&read_cif &get_1st_model);

use strict;
use Storable qw/ nstore retrieve/;


my $cifcachedir="/tmp";

$cifcachedir = "$ENV{HOME}/data/.cif.cache" if defined $ENV{HOME}; 

# hash structure:
# h->{name} = "1je8" (pdb / data_block name)
# |->{pdbchain}->{$pdbchain} = $entryid
# `->{chain}
#       |
#       `->{$entryid}
#             |
#             |->{type}
#             |->{seq}
#             |->{site}->{$act_site_name}->[seqpos1, seqpos2, .., seqposN]
#             |
#             |->{pdbmap}
#             |     |->{ch}->{$seqpos}->[pdbchain1, pdbchain2, ...]
#             |     `->{num}->{$seqpos} = "ALA233"
#             `->{model}
#             |     |
#             |     `->{$pdbmodelid}->{$seqpos}->[..,$i, ... ]
#             |       ^-----_______                  |->{sym} = "O"
#             |                    \                 `->{x,y,z}
#             `->{pdbchain_model}->{A}
#                               |
#                               `->{B}->{$pdbmodelid}->{$seqpos}->[..{x y z }..]

sub get_1st_model{
    my ($cif, $ch) = @_;
    my $h = $cif->{chain}{$ch}{model};
    return undef unless defined $h;

    $h->{(sort {$a <=> $b} keys %$h)[0]};
}

sub read_cif{
    my $file = shift;
    my %p = @_;
    my $all_models = $p{allmodels};
    my $get_atoms  = !$p{noatoms};
    my $cachedir   = defined $p{cachedir} ? $p{cachedir} : $cifcachedir;
    my $cachefile;
    $cachefile = "$cachedir/$1" if $file =~ /([^\/]+)$/;
    if(!$p{cacheoff} && -e $cachefile){
	return retrieve($cachefile);
    }

    my $fh;
    if(ref($file) eq 'GLOB'){ 
	$fh = $file; 
	$p{cacheoff} = 1;
    }
    else { 
	$file = "gunzip < $file |" if $file =~ /\.(?:Z|gz)$/ 
	    && $file =~ /^([-\/\.\w\@ ]+)$/;
	    
	open($fh, $file) || die "Can not open file: $!";
    }

    my %h;
    my $struct;

    local $/ = "\ndata_"; 
    my $cif = <$fh>;

    $h{name} = lc $1 if $cif =~ /^(?:data_)?(\S+)/;

    my $get_chseq = sub{
	my ($t,$r,$ch) = @_; 
	$h{chain}{$ch}{type} = $r->{"$t.type"}; 
	$h{chain}{$ch}{seq}  = $r->{"$t.pdbx_seq_one_letter_code_can"}; 
	$h{chain}{$ch}{seq}  =~ s/\s+//g;
    };

    # get source sequences
    process_std_table($cif, "_entity_poly", $get_chseq);

    # XXX check is obsolete?
    unless( exists $h{chain} ){
	my $r = read_vars(\$cif, "^_entity_poly\\.");
	$get_chseq->("_entity_poly", $r, $r->{"_entity_poly.entity_id"});
    }


    # get pdb mappings
    my %rmap; # reverse pdbmap
    process_std_table($cif, "_pdbx_poly_seq_scheme", 
	sub{
	    my ($t,$r,$ch) = @_; 
	    if( $r->{"$t.pdb_seq_num"} !~/[^\d\-]/ 
		# && $r->{"$t.pdb_seq_num"} >= 0  # they can be < 0 (8-0)
		&& $r->{"$t.pdb_mon_id"} !~ /^n\/a$/i){
		    my $pos = $r->{"$t.seq_id"};
		    my $chn = $r->{"$t.pdb_strand_id"};
		    my $p   = \%{$h{chain}{$ch}{pdbmap}};

		    my  $res;
		    $p->{num}{$pos} = $res = 
			join("",@{$r}{"$t.pdb_mon_id","$t.pdb_seq_num"});

		    $rmap{$ch}{$res} = $pos;

                    unless($chn eq '?' || $chn eq '.'){
			push @{$p->{ch}{$pos}}, $chn;
			$h{pdbchain}{$chn} = $ch;
		    }

	    }
	}
    );

    
    # get atom coordinates
    if($get_atoms){
	my %model = ();

	process_std_table($cif, "_atom_site", 
	sub{
	    my ($t,$r) = @_; 
	    my $ch = $r->{"$t.label_entity_id"};
	    my $m  = $r->{"$t.pdbx_PDB_model_num"};
	    my $pch= $r->{"$t.auth_asym_id"};

	    unless($all_models){
		defined $model{$ch} &&  $m!=$model{$ch} && return;
		$model{$ch} = $m;
	    }


	    my $x = \@{$h{chain}{$ch}{pdbchain_model}{$pch}{$m}{$r->{"$t.label_seq_id"}}};

	    push @$x,
	    {
		(map { $_ => $r->{"$t.Cartn_$_"} } qw(x y z)),
		sym => $r->{"$t.type_symbol"},
	    };

	}
	);

	foreach my $ch (keys %{$h{chain}}){
	    my $t = $h{chain}{$ch}{pdbchain_model};
	    $h{chain}{$ch}{model} = $t->{(sort keys %$t)[0]};
	}
    }

    # get active site
    my $d = 1;
    my %sites = (); # flags duplicated SITEs
    process_std_table($cif,"_struct_site_gen",
	sub{
	    my ($t, $r) = @_;
	    my $chid=$h{pdbchain}{$r->{"$t.auth_asym_id"}};
	    defined $chid || return;
	    my $ch = $h{chain}{$chid};
	    my $res = $r->{"$t.auth_comp_id"}.$r->{"$t.auth_seq_id"};
	    if(defined $ch){
		push @{$ch->{site}{$r->{"$t.site_id"}}},
		    $rmap{$chid}{$res} unless $sites{$rmap{$chid}{$res}}++;
	    }
	}
    );

    # get pdb helix

#     process_std_table($cif,"_struct_conf",
# 	sub{
# 	}
#     );

    # process_std_table($cif,"_struct_conf.pdbx_PDB_helix_length",
    # sub{
    # };

    nstore \%h, $cachefile unless $p{cacheoff};
    return \%h;
}

sub process_std_table{
    my ($t,$sub) = @_[1,2];
    if( my $next_rec = get_loop($_[0], "^$t\\.") ){
	while(my $r = $next_rec->()){
	    my $ch = $r->{"$t.entity_id"};
	    $sub->($t,$r,$ch);
	}
    }
}

sub start_read_loop{
    my $sref = \$_[0];

    
    return sub{
	while($$sref =~/^loop_\s*((?:^_.+?\n)+[^_][\S\s]+?)(?:(?=^(?:_|loop_))|\Z)/mg){
	    my $s = $1;
	    next if defined $_[0] && !($s=~/$_[0]/);

	    return start_read_rec($s);
	}
	undef;
    }
}

sub get_loop{
    return undef unless defined $_[0] && defined $_[1]; 

    pos($_[0]) = 0; # search from the beginning
    # this regexp is crafted specially to avoid perl v5.8.4 segfault
    # XXX do not simplify 
    while($_[0] =~ /^loop_\s*((?:^_.+?\n)+[^_][\S\s]+?)(?:(?=^(?:_|loop_))|\Z)/mg){
	my $s = $1; 
	if($s =~/$_[1]/){ 
	    return start_read_rec($s);
	}
    }
    if($_[0]=~/$_[1]/m){
	my $rs = \$_[0];
	my $re =  $_[1]; 
	my $cnt = 0;
	return sub{ $cnt++ || return read_vars($rs, $re); undef};
    }
    undef;
}

sub read_values{
    if($_[0]=~ /^;((?:.*\n)+?)^;|^([^#_;].+\n)/mg){
	#have to take care of special case of quotingx like: ATOM 'O5'' ..
	return (defined $1 ? $1 : grep {defined && /./}
		split /['"](.*?)['"](?:\s+|\Z)|\s+/, $2)
    }
    ();
}

sub start_read_rec{
    my $s = shift;
    my @fields; 
    my @x = ();
    push @fields,$1 while $s =~ /\G(_\S+)\s+/g;

    my $r = sub {
	while(@x<@fields){ 
	    my @t=read_values($s);
	    @t>0 && push(@x,@t) || last
	} 
	@x>=@fields && return {map {$fields[$_] => shift @x} (0..$#fields)};
    };
}

sub read_vars{
    my ($rs, $re) = @_;
    $re = "" unless defined $re;
    my %h;
    while( $$rs=~/^(_.+?)^(?=(?:_|loop_)|\Z)/mgs ){
	my $s = $1;
	next unless $s=~/$re/;
	$s =~s/^(_\S+)\s+//;
	my $var = $1;
	$h{$var} = (read_values($s))[0];
    }
    \%h;
}

1;
