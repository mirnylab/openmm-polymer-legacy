#!/usr/bin/perl -w
# $Id: knotsubmit.pl,v 1.26 2007/07/23 20:32:46 kolesov Exp $
use strict;
use CGI qw/ :standard :cgi-lib/;
use File::Temp qw/ tempfile tempdir /;
# use IPC::Open2;
use CGI::Carp qw(fatalsToBrowser set_message);
use Storable qw/ nstore retrieve/;
BEGIN {
    sub handle_errors {
	my $msg = shift;
	print "<h1>Error</h1>";
	print "<p>Script returned error: $msg</p>";
    }
    set_message(\&handle_errors);
}


my $CIFPATH = '/net/bio/db/rcsb/all/mmCIF';
my $PDBPATH = '/net/bio/db/rcsb/all/pdb';
my $CGIPATH = '/var/www/cgi-bin/knots';
my $TMP     = "/tmp/cgiknot-".getpwuid($<);
my $HROOT   = "/knots";
my $pdb2cif = "$CGIPATH/pdb2cif.pl";
my $analyzestructure="$CGIPATH/analyze_structure_kgb.pl";
my $pdbref = "getcif.pl";
my $redref = "getreduced.pl";
my $flank  = 8;

my $mystring = query_string();
$mystring="" unless defined $mystring;

umask 022;
unless(-d $TMP){
    mkdir $TMP || die "Can not create temp dir: $!";
}

if(defined param("results")){
    my $fname = sanitize_meta(param("results"));
    #$fname = "$TMP/$fname" if $fname !~ /^\//;
    check_file_indir($TMP,$fname) || die "wrong file name supplied or file expired";
    my $h = retrieve($fname);
    # print h2($fname);
    print "Content-Type: application/x-zip\n\n";

    pack_result($h->{result}, $h->{ciffile}, $h->{cif_base}, $h->{hpdbid});
    exit(0);
}

my %nknot = ( 
    "3_1" => {n => '3 1', name => 'trefoil',      img => "$HROOT/3_1.jpg"},
    "4_1" => {n => '4 1', name => 'figure eight', img => "$HROOT/4_1.jpg"},
    "5_1" => {n => '5 1', name => undef,          img => "$HROOT/5_1.jpg"},
    "5_2" => {n => '5 2', name => undef,          img => "$HROOT/5_2.jpg"},
    "xxx" => {n => '? ?', name => undef,          img => "$HROOT/question.png"},
);

map {
    $_->{un} = $_->{n}; $_->{un} =~ s/ /_/; $_->{n} =~ s/ (\d+)$/<sub>$1<\/sub>/
} values %nknot;


my $email = "grigory \/at\\  mit [dot] edu";
my $logoleft="/knots/knot_full_small.jpg";
my $logoright="/knots/knot_small.jpg";
my %Style = (-title => "Knots in the proteins - prediction server", 
    -author=>$email, -BGCOLOR=>'white', -TEXT=>'black');

print header;
print start_html(%Style);


if(defined param("jmol_results")){
    my ($fname,$i) = param("jmol_results");
    $fname = sanitize_meta($fname);
    $fname = "$TMP/$fname" if $fname !~ /^\//;
    check_file_indir($TMP,$fname) || die "wrong file name supplied or file expired";

    my $h = retrieve($fname);
    $i>=0 && $i<scalar(@{$h->{result}}) && scalar(@{$h->{result}})>0 || die "wrong counter or file name supplied";

    start_script_head(scalar(@{$h->{result}}), $i);
    show_jmol($h->{result}[$i], $h->{ciffile}, $h->{cif_base}, $i,$h->{hpdbid});
    print end_html; exit(0);
}
elsif(defined param("pdbid") || defined param("uploaded_file")){
    #start_script();

    my ($cif_fh, $cif_n) = tempfile(DIR=>$TMP, SUFFIX=> ".cif");
    my $pdb_n;
    my $cif_base = $cif_n;
    my $iscif = "";
     

    if(defined param("pdbid") && param("pdbid") !~ /^\s*$/){
	my $pdbid = sanitize_meta(param("pdbid"));
	close $cif_fh;
	#system("gunzip < $CIFPATH/$pdbid.cif.Z > $cif_n")==0 
	system("gunzip < $PDBPATH/pdb$pdbid.ent.Z > $cif_n")==0 
	    || do{
		print "<h3>PDB entry $pdbid not found.</h3>";
		print end_html(); exit(0);
	    }; 
	    # $iscif = "-c";
    }
    else{
	my $up = upload("uploaded_file");
	defined $up || do{
		print "<h3>Error: please enter PDB id or upload a structure.</h3>";
		print end_html(); exit(0);
	};
	my $first = <$up>;
	print $cif_fh $first;
	while(<$up>){ print $cif_fh $_};
	close $up;
	close $cif_fh;

	$iscif = "-c" unless $first =~ /^(?:ATOM|HETATM|REMARK|HEADER)/;
    }


    my $straight = "";
    $straight = "-s" if param('straight') eq 'on';
    open(my $rh, "$analyzestructure $straight $iscif -r $cif_n $cif_n |");
    # open(my $rh, "$analyzestructure $cif_n |");

    # open(my $errh, ">$TMP/knot.err");
    my $out = join '',<$rh>;
    close $rh;

#     if(defined $pdb_n){
# 	unlink $cif_n;
# 	$cif_n = $pdb_n;
#     }

    display_result($out, $cif_n, $cif_base);
    #system("/bin/echo \"/bin/rm -f $cif_n\" | /usr/bin/at \"now + 10 minutes\" > /dev/null 2>&1"); 

    print end_html; exit(0);
}


# home default
print "<center>";
print  img({-src=> $logoleft,-align=>'absmiddle'}),
"<font size=50>Protein Knots</font>",
img({-src=> $logoright, align=>'absmiddle'});
print "</center>";


print <<ABS;
The knot server allows the user to check pdb entries or uploaded structures for
knots and to visualize them. The size of a knot is determined by deleting amino
acids from both ends. This procedure is, however, not perfect and the resulting
size should only be treated as a guideline.<p>

ABS

print "<br><p>\n";

print start_multipart_form(-action => "/cgi-bin/knots/knotsubmit.pl"),"\n";

print "Please enter pdb id (e.g. 1v2x): ",br, textfield('pdbid'),p,"\n";
print "Or upload file (pdb or mmCIF):",br,"\n";
print filefield(-name=>'uploaded_file',
    -default=>'starting value',
    -size=>50,
    -maxlength=>80),"\n",br;

# print checkbox(-name=>'straight',
#     -checked=>1,
#     -value=>'ON',
#     -label=>'Connect gaps in structure by a straight line'),p;


print radio_group(-name=>'straight',
    -values=>['on','off'],
    -default=>'on',
    -linebreak=>'true',
    -labels=>{
	on => "Connect gaps in the structure by straight line",
	off=> "Treat each unconnected fragment separately",
    },
),p;





print submit(-value => "Submit", -name => "submit"),p,"\n";
print end_form(),"\n";

if(param("showlist")){
    my $t = $mystring;
    $t =~ s/showlist[^&;+]+[&+;]?//;
    print "<a href='/cgi-bin/knots/knotsubmit.pl?$t'>List of known knots</a><br>";
    show_list();
}
else{
    print "<a href='/cgi-bin/knots/knotsubmit.pl?showlist=1&$mystring'>List of known knots</a><br>";
}


if(param("knot_definition")){
    my $t = $mystring;
    $t =~ s/knot_definition[^&;+]+[&+;]?//;
    print "<a id=def href=/cgi-bin/knots/knotsubmit.pl?$t#def>How we define knots</a>";
    print p,knot_definition();
}
else{
    print "<a id=def href=/cgi-bin/knots/knotsubmit.pl?knot_definition=1&$mystring#def>How we define knots</a>";
}



print <<REF;
<h4>References</h4>
<ul>
<li><a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Retrieve&dopt=AbstractPlus&list_uids=17517776&query_hl=1&itool=pubmed_docsum'>Kolesov G, Virnau P, Kardar M, Mirny LA. Protein knot server: detection of knots in protein structures. Nucleic Acids Res. 2007 May 21</a> <a href='http://nar.oxfordjournals.org/cgi/reprint/gkm312?ijkey=GFbISFYeK059SYi&keytype=ref'>[PDF]</a><br>
<li><a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Retrieve&dopt=AbstractPlus&list_uids=16978047&query_hl=4&itool=pubmed_Brief'>Virnau P, Mirny LA and Kardar M. Intricate knots in proteins: function and evolution. PLoS Comput Biol. 2006 Sep 15;2(9):e122. Epub 2006 Jul 28.</a><br>
<li><a href='http://www.mit.edu/~leonid/index.html'>Mirny lab webpage</a>
</ul>
<p>
Questions about particular knots? Contact <a href="mailto:virnau [at] uni-mainz [dot] de">Peter Virnau</a>.
<hr>
Contact <a href="mailto:$email">webmaster</a>
REF


print end_html;




sub display_result{
    my ($result, $ciffile,$cif_base) = @_;
    #print STDERR $result;
    # print pre($result);

    my @r =  map { $_->{chn} =~ s/(?:^'|'$)//g; $_} map { 
	{
	    chn   => $_->[0],
	    knot  => $nknot{$_->[4]},
	    start => $_->[5], 
	    stop  => $_->[6],
	    chn1  => $_->[1],
	    chn2  => $_->[2],
	}
    } map {
	$_->[4]='xxx' unless defined($nknot{$_->[4]}) ||  $_->[4] eq '0_0';
	$_;
    } map {[split]} grep(/^\S+\s+\-?\d+/,(split /\n/, $result));

    
    @r = grep {defined $_->{knot}} @r;
#     foreach my $knot (@r){
# 	$knot->{start} -= $flank if $knot->{start}>$flank+1;
# 	$knot->{stop}  += $flank;
#     }
    

    for(my $i=0; $i<@r; $i++){
	for(my $j=$i+1; $j<@r; $j++){
	    if($r[$i]{chn} eq $r[$j]{chn} &&  $r[$i]{knot} == $r[$j]{knot}  
		&& max($r[$i]{start},$r[$j]{start})<min($r[$i]{stop},$r[$j]{stop}) ){

		$r[$j]{start} = min($r[$i]{start}, $r[$j]{start});
		$r[$j]{stop}  = max($r[$i]{stop},  $r[$j]{stop});
		$r[$i] = undef;
		last;
	    }
	}
    }

    @r = grep {defined} @r;

    my $hpdbid = sanitize_meta(param("pdbid"));
    $hpdbid = "user-defined" if !(defined $hpdbid) || $hpdbid eq "";
    my $rcsbid = rcsb_link($hpdbid);

    unless(@r>0){
	print "<h4>No knots found in the $rcsbid structure.</h4>";
	return;
    }

    print "<h3>Knots found in the $rcsbid structure:</h3>",br;
    my $i = 0;

    my $pass_n=pass_result_file(\@r, $ciffile, $cif_base, $hpdbid);
    
    print table({-border=>1},
	Tr({-align=>'CENTER'},[
	    th(['Knot residues','Chain start-stop','Knot type','Knot', '']),
	    map{
		my $r = $_;
		td(
		    [
		    "$r->{start}-$r->{stop}$r->{chn}",
		    "$r->{chn1}-$r->{chn2}$r->{chn}",
		    "<a href=http://katlas.math.toronto.edu/wiki/Image:$r->{knot}{un}.gif>$r->{knot}{n}</a>".(
			defined $r->{knot}{name} ? 
			    "<br>($r->{knot}{name} knot)" : "",
			),
		    img({-src => $r->{knot}{img}}), show_jmol_form($pass_n, $i++)
		    ]
		)
	    } @r,
	    ]
	)),br;

    # save_result_form(\@r, $ciffile, $cif_base, $hpdbid);
    save_result_form($pass_n);

#     for(my $i=0; $i<@r; $i++){ 
# 	show_jmol($r[$i],$ciffile,$cif_base,$i,$hpdbid);
#     }
#2 
#     print "<br>" x 20;
}

sub save_result_form{
    my ($pass_n) = @_;

    # my $pass_n = pass_result_file(@_);

    print start_multipart_form(-action => "/cgi-bin/knots/knotsout.zip"),"\n";
    print hidden(-name => 'results', -default => $pass_n);

    print submit(-value => "Download results and rasmol scripts as zip package", -name => "save_res"),p,"\n";
    print end_form(),"\n";
}

sub show_jmol_form{
    my ($pass_n,$i) = @_;
    
    my $fbase="\n"; 
    $fbase = $1 if $pass_n =~ /\/([^\/]+)$/;


    my $s="<a href=/cgi-bin/knots/knotsubmit.pl?jmol_results=$fbase&jmol_results=$i>Jmol visualization</a>";
#     $s = start_multipart_form(-action => "/cgi-bin/knots/knotsubmit.pl")."\n";
#     $s.= hidden(-name => 'jmol_results', -default => [$pass_n,$i]);
# 
#     $s.= submit(-value => "Show Jmol visualization", -name => "show_jmol_$i").p."\n";
#     $s.= end_form()."\n";
    return $s;
}

sub pass_result_file{
    my ($r, $ciffile, $cif_base, $hpdbid) = @_;

    my ($pass_fh, $pass_n) = tempfile(DIR=>$TMP, SUFFIX=> ".pas");
    close $pass_fh;
    nstore {
	result => \@$r, 
	ciffile=> $ciffile, 
	cif_base=> $cif_base, 
	hpdbid => $hpdbid
    }, $pass_n;

    return $pass_n;
}


sub pack_result{
    my ($r, $ciffile, $cif_base, $hpdbid) = @_;

    my $tdir = tempdir(DIR=>$TMP, CLEANUP=>1);
    my $sname = $hpdbid;
    $sname =~ s/ /_/;

    my $dirname = "$tdir/Knots_in_$hpdbid";
    mkdir $dirname || die "Can not create tempdir : $!";

    open(my $fr, "<", $ciffile)       || die "Can not open file: $!";
    open(my $fw, ">", "$dirname/$sname") || die "Can not open file: $!";
    while(<$fr>){print $fw $_};
    close $fr;
    close $fw;

    for(my $i=0; $i<@$r; $i++){ 
	create_scripts($dirname,$r->[$i],$sname,$cif_base,$i);
    }

    save_readme("$dirname/README");
    save_results_table($r,"$dirname/result_knots.txt");

    chdir $tdir || die "Can not chdir : $!";
    system("zip -qr - Knots_in_$hpdbid"); 

}


sub create_scripts{
    my ($dirname,$r,$ciffile,$cif_base,$i) = @_;

    open(my $fh, ">", "$dirname/knot_$i.spt");
    print $fh "load \"$ciffile\"\n";
    print $fh join("\n",(ras_script($r)));
    close $fh;

#     open(my $fr, "<", "$cif_base.$r->{chn}-$r->{chn1}") 
# 	|| die "Can not open file: $!";

    open(my $fw, ">", "$dirname/knot_$i"."_red.pdb") 
    || die "Can not open file in $dirname: $!";
    # while(<$fr>){print $fw $_};
    getreduced_pdb("$cif_base.$r->{chn}-$r->{chn1}", $fw);
#   close $fr;
    close $fw;

    

    open ($fh, ">", "$dirname/knot_$i"."_red.spt") || die "Can not create file: $!";
    print $fh "load \"knot_${i}_red.pdb\"\n";
    print $fh join("\n",(grep !/zoom/,red_script(100)));
    close $fh;

}

sub getreduced_pdb{
    my ($file, $outfh) = @_;
    open(my $fh, "<", "$file") || die "Can not open file: $!";
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

	printf $outfh "%-6s %4d %1s%-3s %3s %1s %3d    %8.3f%8.3f%8.3f  %4.2f%6.2f %11s\n", 'ATOM', $i, '', 'CA', 'GLY', 'A', $i, $t->[1], $t->[2], $t->[3], 1.00, 99.00, 'C';
    }
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

sub rcsb_link{
    my $id = shift;
    return 
    "<a href=http://www.rcsb.org/pdb/explore/explore.do?structureId=".uc($id).">$id</a>"
	if length $id == 4;
    $id;
}

sub show_jmol{
    my ($r, $ciffile, $cif_base, $ix,$hpdbid) = @_;
    my $cifname = $ciffile;
    $cifname =~ s/\Q$TMP\/\E//;
    $cif_base =~ s/\Q$TMP\/\E//;


    # my $scrprelude = "cpk off; wireframe off; cartoon; color gray;"; 

    my ($lstart, $lstop) = ($r->{start} - $flank, $r->{start}-1);
    my ($rstart, $rstop) = ($r->{stop}  + 1, $r->{stop}+$flank);
    $lstart = 1 if $lstart<1;

    my $rasscript = ras_script($r);
    my $redscript = red_script();

    print "<h3><a id=$ix>Residues $r->{start}-$r->{stop}$r->{chn}</a></h3><br>\n";
#     print 'Type of knot found:',br;
#     print table({-border=>1},
# 	Tr({-align=>'CENTER'},[
# 	    td(["<a href=http://katlas.math.toronto.edu/wiki/Image:$r->{knot}{un}.gif>$r->{knot}{n}</a>".(
# 	defined $r->{knot}{name} ? "<br>($r->{knot}{name} knot)" : "",
# 		), img({-src => $r->{knot}{img}}) ]
# 	    )
# 	    ]
# 	)),br;

    my $rcsbid = rcsb_link($hpdbid);


print <<JAPPLET;

<table border=1>
<tr>
   <th>Knot in the $rcsbid structure</th><th>Simplified representation of the knot</th>
</tr>
<tr>
   <td width=350 align=left>
 
    <applet name="jmol2_$ix" code="JmolApplet" archive="/knots/JmolApplet.jar"
	width="350" height="350" align="left" mayscript="true">
	<param name="load"    value="$pdbref?$cifname">
	<param name='progressbar' value='true'>
	<param name="bgcolor" value="#ffffff">
	<param name="style"   value="wireframe">
	<param name="width"   value="350">
	<param name="height"  value="350">
	<param name="wireframeRotation" value="true">
	<param name="perspectiveDepth"  value="true">
	<param name="script" value="$rasscript;">
    </applet>
</td>
<td width=350 align=left>
    <applet name="jmol2_r$ix" code="JmolApplet" archive="/knots/JmolApplet.jar"
	width="350" height="350" align="left" mayscript="true">
	<param name="load"    value="$redref?$cif_base.$r->{chn}-$r->{chn1}">
	<param name='progressbar' value='true'>
	<param name="bgcolor" value="#ffffff">
	<param name="style"   value="wireframe">
	<param name="width"   value="350">
	<param name="height"  value="350">
	<param name="wireframeRotation" value="true">
	<param name="perspectiveDepth"  value="true">
	<param name="script" value="$redscript;">
    </applet>
</td>

</tr></table>

<table border=0>
<tr><td>
<form name="jtoggle_$ix" onSubmit="
 toggle_jmol2_$ix('($r->{start}-$r->{stop}$r->{chn}) or ($lstart-$lstop$r->{chn}) or ($rstart-$rstop$r->{chn})');return false;">
<input type="submit" name=btoggle_$ix Value="Hide/Show unknotted structure">
</form>
</td>

<td>
<form name="jspin_$ix" onSubmit="spin_jmol2_$ix();return false;">
<input type="submit" name=bspin_$ix Value="Spin the structures">
</form>
</td></tr></table>

<pre>
Hint: hold Ctrl-Alt to move the structure, Shift to zoom. Right click to get console.
Knotted region is defined as 'knot', typing 'select knot' will select corresponding residues.
</pre>
Enter one-line RasMol/Chime script commands here:
<form name="line_$ix" onSubmit="
 run_jmol2_$ix(document.line_$ix.script_$ix.value);document.line_$ix.script_$ix.select();return false;">
<input name="script_$ix" onFocus="document.line_$ix.script_$ix.select()" size=80 type=text>
</form>

<hr>
<p>

JAPPLET

}

sub ras_script{
    my ($r) = @_;

    my ($lstart, $lstop) = ($r->{start} - $flank, $r->{start}-1);
    my ($rstart, $rstop) = ($r->{stop}  + 1, $r->{stop}+$flank);
    $lstart = 1 if $lstart<1;

    my @ret =  (
       "cpk off", "wireframe off","cartoon","color gray",
       "select ligand,rna,dna",
       "color chain",
       "wireframe 100", 
       "define knot ($r->{start}-$r->{stop}$r->{chn}) or ($lstart-$lstop$r->{chn}) or ($rstart-$rstop$r->{chn})", 
       "select knot",
       "color group", 
       "cartoon off", 
       "backbone 200",
       "center selected"
    );

    return @ret if wantarray;
    return join(";",@ret);
}

sub red_script{
    my ($backbone_d) = @_;
    $backbone_d = 5 unless defined $backbone_d;
    my @ret = ("wireframe off", "cpk off", "backbone $backbone_d","color group","zoom 400");

    return @ret if wantarray;
    return join(";", @ret);
}

sub start_script_head{
    my ($n,$ix) = @_;
    my $once = -1;
    if(defined $ix){
	$once = 1;
    }
    else{ $ix = 0;}

    print "<script>\n";
    for(my $i=$ix; $i<$n && $once--; $i++){
    print <<FUNC;
function run_jmol2_$i(script) {
  document.jmol2_$i.script(script);
}
jmol2_${i}_spin = false;
function spin_jmol2_$i(script) {
  if(jmol2_${i}_spin){
      jmol2_${i}_spin = false;
      document.jmol2_$i.script("spin off");
      document.jmol2_r$i.script("spin off");
  }
  else{
      jmol2_${i}_spin = true;
      document.jmol2_$i.script("spin 20");
      document.jmol2_r$i.script("spin 20");
  }
}
jmol2_${i}_s_hide = false;
function toggle_jmol2_$i(selectres) {
  if(jmol2_${i}_s_hide){
     jmol2_${i}_s_hide = false;
     document.jmol2_$i.script("display all;");
  }
  else{
     jmol2_${i}_s_hide = true;
     document.jmol2_$i.script("select "+selectres+"; display selected;");
  }
}
FUNC
    }

    print q|
function br() {
 document.write("<br>");
}
</script>
|;


}

sub show_list{
    my @list = ();
    # my $straight = "";
    # $straight="&straight=on" if param("straight");
    while(<DATA>){
	my @t = map {s/^\"(.*?)\s*\"$/$1/;$_} split /\t+/;
	$t[2] = "<a href=/cgi-bin/knots/knotsubmit.pl?pdbid=$t[2]&straight=on>$t[2]</a>" if defined $t[2];
	$t[4] =~ s/(\d+) (\d+)/$1<sub>$2<\/sub>/;
	push @list, \@t;
    }

    # print "List of known knots:<br>";

    print table({-border=>0},
	Tr({-align => "LEFT", -valign=>"CENTER"},
	    [
	    th(['Protein','Species','PDB code','Length','Knot','Knotted core']),
	    ( map {td($_)} @list )
	    ]
	)
    );
    close DATA;
}

sub knot_definition{
    return 
'

Mathematically, <a href=http://mathworld.wolfram.com/Knot.html>knots</a> are
only well defined in closed (circular) loops. However, both the N- and
C-termini of open proteins are typically located close to the surface of the
protein and can be connected unambiguously: We reduce the protein to its
backbone and draw two lines outward starting at the termini in the direction of
the connection line between the centre of mass of the backbone and the
respective ends. The two lines are joined by a big loop, and the structure is
topologically classified by the determination of its Alexander polynomial. To
determine an estimate for the size of the knotted core, we successively delete
amino acids from the N-terminus until the protein becomes unknotted. The
procedure is repeated at the C-terminus starting with the last N-terminal
deletion structure that contained the original knot. For each deletion, the
outward-pointing line through the new termini is parallel to the respective
lines computed for the full structure. Unfortunately, the size of a knot is not
always precisely determined by this procedure, so reported sizes should
therefore only be treated as approximate.<p>

To speed up calculations, the KMT reduction scheme is used. This algorithm
successively deletes amino acids that are not essential to the topological
structure of the protein. It is also employed to create a reduced
representation of the knot.  In the course of our investigations we came up
with a number of stringent criteria that a structure should satisfy to be
classified as knotted:<p>

<ol>
<li>The Alexander polynomial should yield a knot. </li>
<li>There should not be any gaps in the polypeptide backbone. (See below.) </li>
<li>The knot should persist if two amino acids are removed from each end. (This prevents knots formed by just a few residues at the end of the chain passing through the loop - "shallow knots" and knots which only appear due to our specific loop closure procedure.)</li>
</ol>
<p>

Unfortunately, there are some structures containing regions of the backbone
that were not resolved and for which coordinates are not reported in PDB (a gap
in the structure). Mobile loops may not be resolved by X-ray crystallography
unless they are stabilized by a ligand or by protein engineering, for example.
If the polypeptide chain contains a gap, the knot is reported if a knot is
present in at least one fragment of the chain and  the structure that
results from gaps being bridged with straight lines contains a knot. These
criteria form the basis of our list of known knots. We have also included
knotted structures with gaps if at least one homolog is knotted.
'

}

sub save_readme{
    my ($file) = @_;
    open(my $fh, ">", $file) || die "Can not open file: $!";
    print $fh <<README;
Rasmol scripts found in this directory are called knot_X.spt for the knots in
the structure and knot_X_red.spt for reduced version of the knots (where X is
number of the knot starting with 0).  You can run them by typing in the command
line:

rasmol -script knot_d_0.spt

The summary of results is found in file called "result_knots.txt";

Rasmol does not understand cif format. For cif files you should use Chime, Jmol
or similar.

README

    close $fh;
}

sub save_results_table{
    my ($res, $file) = @_;
    open(my $fh, ">", $file) || die "Can not open file: $!";
    my $oldsel = select $fh;
    my $head = sprintf "%5s %s-%s %9s\n", "Chain","start","stop","Knot type";
    print $fh $head;
    print $fh "-" x (length($head)-1),"\n";
    foreach my $r (@$res){
	printf $fh "%-5s %5d-%-4d %-9s\n", $r->{chn}, $r->{start}, $r->{stop}, $r->{knot}{un};
    }
    close $fh;
    select $oldsel;
}


sub reformat_line{
    my ($line,$atomn) = @_;
    my @x = split /\s+/,$line;
    my $s=join " ",('ATOM', ++$$atomn, $x[1], $x[2],'.', $x[3], $x[4],1,$x[5],'?',@x[7..9]);
    $s.="\n";
    return $s;
}



sub sanitize_meta{
    $_[0] =~ s/(\\+)([\$\@`'"\t\n\r\^*;!><|\{\}\[\]\(\)~?& ])/('\\' x (length($1) + (length $1)%2))."$2"/eg;
    $_[0] =~      s/([\$\@`'"\t\n\r\^*;!><|\{\}\[\]\(\)~?& ])/\\$1/g;
    $_[0];
}

sub check_file_indir{
    my ($dir, $file) = @_;

    my $fbase="\n"; 
    $fbase = $1 if $file =~ /\/([^\/]+)$/;
    opendir(my $dh, $dir) || die "Can not read dir: $!";
    my $ret  = 0;

    if($file !~ /\/\.\.\// && grep /^$fbase$/,readdir($dh) ){ $ret = 1}
    close $dh;

    return $ret;
}


sub min{
  my $min = $_[0];
  foreach my $k (@_){   $min=$k if defined($k) && (!defined($min) || $k<$min);  }
  return $min;
}

sub max{
  my $max = $_[0];
  foreach my $k (@_){ $max=$k if defined($k) && $k>$max;  }
  return $max;
}
# ""	"Mus musculus "	"2znc "	249	"3 1 "	"32-246 (3) 1 "

__DATA__
"YbeA-like "	"E.coli "	"1ns5 "	153	"3 1"	"69-121 (32) "
""	"T.maritima "	"1o6d "	147	"3 1 "	"68-117 (30) "
""	"S.aureus "	"1vh0 "	157	"3 1 "	"73-126 (31) "
""	"B.subtilis "	"1to0 "	148	"3 1 "	"73-125 (32) "
"tRNA(m1G37)-methyltransferase TrmD "	"H.influenza "	"1uaj "	241	"3 1 "	"85-130 (92)"
""	"E.coli "	"1p9p "	235	"3 1 "	"90-130 (89)"
""	"S.cerevisiae"	"2v3k "	219	"3 1 "	"175-225 (27)"
"SpoU-like RNA 2'-O ribose mtf. "	"T.thermophilus "	"1v2x "	191	"3 1 "	"96-140 (51) "
""	"H.influenza "	"1j85 "	156	"3 1 "	"77-114 (42) "
""	"T.thermophilus "	"1ipa "	258	"3 1 "	"190-234 (29)"
""	"E.coli "	"1gz0 "	242	"3 1 "	"173-215 (28) "
""	"A. aeolicus "	"1zjr "	197	"3 1 "	"100-144 (58) "
""	"S. viridochromog. "	"1x7p "	267	"3 1 "	"209-251 (31)"
"YggJ C-terminal domain-like "	"H.influenza "	"1nxz "	246	"3 1 "	"166-217 (30) "
""	"B.subtilis "	"1vhk "	235	"3 1 "	"168-226 (27) 1 "
""	"T.thermophilus "	"1v6z "	227	"3 1 "	"104-203 (25) 3 "
"Hypothetical protein MTH1 (MT0001) "	"A.M. thermoautotr. "	"1k3r "	262	"3 1 "	"48-234 (28) "
""	""	""	""	""	""
"Carbonic Anhydrases: "	""	""	""	""	""
"Carbonic Anhydrase "	"N. gonorrhoeae "	"1kop "	223	"3 1 "	"39-226 (0) "
"Carbonic Anhydrase I "	"H.sapiens "	"1hcb "	258	"3 1 "	"29-256 (2) "
"Carbonic Anhydrase II "	"H.sapiens "	"1lug "	259	"3 1 "	"31-257 (3) "
""	"Bos Taurus "	"1v9e "	259	"3 1 "	"32-256 (3) "
""	"Dunaliella salina "	"1y7w "	274	"3 1 "	"39-272 (4) "
"Carbonic Anhydrase III "	"Rattus norv. "	"1flj "	259	"3 1 "	"31-258 (3) "
""	"H. sapiens "	"1z93 "	263	"3 1 "	"31-258 (9) "
"Carbonic Anhydrase IV "	"H.sapiens "	"1znc "	262	"3 1 "	"31-258 (1) "
"Carbonic Anhydrase V "	"Mus musculus "	"1keq "	238	"3 1 "	"31-258 (4) "
"Carbonic Anhydrase VII "	"H.sapiens "	"1jd0 "	260	"3 1 "	"31-258 (3) "
"Carbonic Anhydrase XIV "	"Mus Musculus "	"1rj6 "	259	"3 1 "	"31-258 (2) "
""	""	""	""	""	""
"Miscellaneous: "	""	""	""	""	""
"Ubiquitin Hydrolase UCH-L3 "	"H.sapiens "	"1xd3 "	229	"5 2 "	"13-212 (12) 4 "
""	"S.cerevisiae (synth.) "	"1cmx "	214	"3 1 "	"14-228 (6) 4,1 "
"Ubiquitin Hydrolase UCH-L1 "	"H.sapiens "	"2etl "	219	"5 2 "	"10-216 (13)"
"S-adenosylmethionine synthetase "	"E.coli "	"1fug "	383	"3 1 "	"33-260 (32) "
""	"rattus norv. "	"1qm4 "	368	"3 1 "	"46-281 (29) 1 "
""	"H. sapiens"	"2p02"	380	"3 1"	"59-302 (21)"
"Class II ketol-acid reductoisomerase "	"Spinacia oleracea "	"1yve "	513	"4 1 "	"321-533 (62) "
""	"E.coli "	"1yrl "	487	"4 1 "	"222-437 (52) 2 "
"Transcarbamylase "	"B.fragilis "	"1js1 "	324	"3 1 "	"169-267 (57) "
""	"X.campestris "	"1yh0 "	328	"3 1 "	"173-277 (57) 2 "
"Methyltransferase"	"H.sapiens"	"2ha8"	159	"3 1"	"103-148 (30)"
""	"P. gingivalis"	"2i6d"	231	"3 1"	"177-223 (9)"
