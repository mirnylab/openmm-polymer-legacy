#!/usr/bin/perl
use CGI;
use File::Temp qw/ tempfile/;

use CGI::Carp qw(fatalsToBrowser set_message);
BEGIN {
  sub handle_errors {
      my $msg = shift;
      print "<h1>Error</h1>";
      print "<p>Script returned error: $msg</p>";
  }
  set_message(\&handle_errors);
}

my $CGIPATH = '/var/www/cgi-bin/knots';
my $CIFPATH = '/net/bio/db/rcsb/all/mmCIF';
my $TMP = "/tmp";
my $ANALYZER = "$CGIPATH/analyze_structure.pl";

my $q = new CGI;

if ($q->param()) {
  my $pdbid = sanitize_meta($q->param('pdbid'));
  my $pdbfile = "$CIFPATH/$pdbid.cif.Z";
  print STDERR "$pdbfile\n";
  if ($pdbfile =~ /\.[zZ]$/) {
    open(FD, "gunzip -c $pdbfile |")
      or die "Fatal: can't open input file '$pdbfile'\n";
  } else {
    open(FD, $pdbfile)
      or die "Fatal: can't open input file '$pdbfile'\n";
  }

  

  my ($tfh, $tname) = tempfile(DIR => $TMP, SUFFIX => ".cif");
  while (<FD>) {
    print $tfh $_;
  }
  close $tfh;

  print $q->header(),
        $q->start_html("mmCIF: $pdbid"),
        $q->h1("mmCIF: $pdbid");

  my $results = join '',`$ANALYZER $tname`;
  print STDERR "$ANALYZER $tname\n";
  print $q->h1('Results:'), 
        $q->pre($results),
        $q->end_html();

  unlink $tname;
}

sub sanitize_meta{
    $_[0] =~ s/(\\+)([\$\@`'"\t\n\r\^*;!><|\{\}\[\]\(\)~?& ])/('\\' x (length($1) + (length $1)%2))."$2"/eg;
    $_[0] =~      s/([\$\@`'"\t\n\r\^*;!><|\{\}\[\]\(\)~?& ])/\\$1/g;
    $_[0];
}

