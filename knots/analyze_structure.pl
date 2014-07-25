#! /usr/bin/perl
use strict;
use File::Temp qw/ tempfile tempdir /;
use POSIX;
use Getopt::Long;

# unlink "knot.dat";

my $redknotbase;

GetOptions(
    "r=s" => \$redknotbase,
);

#----------------------------------------------------------------------------------------
# OUTER LOOP: Aquire list of proteins to be analyzed and loop over all mmCIF structures |
#----------------------------------------------------------------------------------------

my $CGIPATH = '/var/www/cgi-bin/knots';
my $CIFPATH = '/net/bio/db/rcsb/all/mmCIF';
#my $TMP = "$CGIPATH/tmp";
my $TMP = "/tmp";

print "START\n";
# @protein_array=`cat complete_list_2006.txt`;
my @protein_array =  ($ARGV[0]);



my @configs = ();

for(my $counter=0 ; $counter < @protein_array; $counter++)
{
  print "LOOP=$counter\n";


  #----------------------------------------
  # 1. Get C-alpha coordinates from mmCIF |
  #----------------------------------------

  #my @coord = `awk '\$1==\"ATOM\" && \$4==\"CA\" {print \$5,\$8,\$9,\$11,\$12,\$13} END {print 0,0,0,0,0,0}' $protein_array[$counter]`;

  my @coord = ();

  {
    while(<>){
      my @x = split;
      push @coord, [@x[4,7,8,10,11,12,6]] if($x[0] eq 'ATOM' && $x[3] eq 'CA');
      last if eof;
    }
    push @coord, [qw/0 0 0 0 0 0/];


  }
  


  # $het := "$5"                 information regarding microheterogeneity (.,A,B,...)
  # $ent := "$8"                 signals when a new structure starts (1,2,...)
  # $ent := "$9"                 sequence number of amino acid
  # $x   := "$11","$12","$13"    x,y,z coordinates

  #------------------------------------------------  
  # 2. Create configuration files for Reduce_knot |
  #------------------------------------------------

  my $monomer_ctr       = 0; #atom number in structure 
  my $configuration_ctr = 0; #configuration file number (k.0,k.1,...)
  my $start             = 0;



  my ($check1, $check2) = (INT_MIN,INT_MIN);
  my $prev_pdbchn = "";
  my (@het, @x, @y, @z);
  my @aa_ix = ();

  foreach my $line (@coord) {
    my ($ent, $aa, $pdbchn);
    ($het[$monomer_ctr],$ent,$aa,$x[$monomer_ctr],$y[$monomer_ctr],$z[$monomer_ctr],$pdbchn) = @$line;
    $aa_ix[$monomer_ctr] = $aa;


    #----------------------------------------------------------------------------------------------
    # 2.1. Split configuration if a new molecule starts or sequence of amino acids is interrupted |
    #      -> print to configuration file                                                         |
    #----------------------------------------------------------------------------------------------

    #check1 and check2 store information of $ent and $aa from last iteration

    if( $monomer_ctr != 0           &&                          #not first line? AND
        ( $ent       != $check1     ||                          #( new molecule? OR
          ($aa       != $check2+1   &&   $aa != $check2) ) ) {  #  interrupted sequence of amino acids) 

      #---------------------------------------------------------------------------------------------
      # 2.1.1. Check for microheterogeneity (file contains multiple coordinates for the same atom) |
      #---------------------------------------------------------------------------------------------

      my $micro_het_symbol_number;
      my @micro_het_symbol;
      my $micro_het_symbol_number = 1;
      $micro_het_symbol[$micro_het_symbol_number] = '.'; 

      for (my $i=1; $i<$monomer_ctr;$i++){
        #symbol already in symbol array ?
	my $k = 0;
        for (my $j=1;$j<=$micro_het_symbol_number;$j++) {
           if($het[$i] eq $micro_het_symbol[$j]) {$k=1;}
        }
        #no: add symbol to symbol array; example: $micro_het_symbol[2]=A
        if($k==0){
           $micro_het_symbol_number++;
           $micro_het_symbol[$micro_het_symbol_number]=$het[$i];
        }
      }

      #print $micro_het_symbol_number; 
      #for ($j=1;$j<=$micro_het_symbol_number;$j++) {
      #   print $micro_het_symbol[$j];
      #}

      #-----------------------------------------------------------------------------------
      # 2.1.2.a. No microheterogeneity: Print coordinates from $start to $monomer_number |
      #-----------------------------------------------------------------------------------

      if($micro_het_symbol_number==1){
   
        #determine number of monomers in configuration file
        my $number_of_monomers=0;
        for (my $i=$start;$i<$monomer_ctr;$i++){
          $number_of_monomers++;
        }

	my ($conf_fh, $conf_file) = tempfile(DIR=>$TMP, SUFFIX=> ".conf");
        print $conf_fh "t=0\n\n$number_of_monomers\n";

        #print coordinates
        for (my ($i,$j)=($start,1);$i<$monomer_ctr;$i++){
          print $conf_fh "$j\t$x[$i]\t$y[$i]\t$z[$i]\n";
          $j++;
        }
   
        close $conf_fh;
	push @configs, {file => $conf_file, chn => $prev_pdbchn,
	  start => $aa_ix[$start]};
        $configuration_ctr++;
      } #end if($micro ...

      #----------------------------------------------------------------------------------------------
      # 2.1.2.b. Microheterogeneity: Print coordinates - each case gets separate configuration file |
      #----------------------------------------------------------------------------------------------

      else{
        for(my $k=2;$k<=$micro_het_symbol_number;$k++) {

          #determine number of monomers in configuration file
          my $number_of_monomers=0; 
          for (my $i=$start;$i<$monomer_ctr;$i++)  {
            if ($het[$i] eq '.' || $het[$i] eq $micro_het_symbol[$k]) {
              $number_of_monomers++;
            }
          }

	  my ($conf_fh, $conf_file) = tempfile(DIR=>$TMP, SUFFIX=> ".conf");

          print $conf_fh "t=0\n\n$number_of_monomers\n";

          #print coordinates
          for (my ($i,$j) =($start,1);$i<$monomer_ctr;$i++)  {
            if ($het[$i] eq '.' || $het[$i] eq $micro_het_symbol[$k]) {
              print $conf_fh "$j\t$x[$i]\t$y[$i]\t$z[$i]\n"; 
              $j++;
            }
          }
	  close $conf_fh;
	  push @configs, {file => $conf_file, chn => $prev_pdbchn,
	    start => $aa_ix[$start]};

	  $configuration_ctr++;
        }
      } 

      $start=$monomer_ctr;

    } #end if( $monomer_ctr != 0 ... -> continue with next protein or protein part in mmCIF

    $monomer_ctr++;
    $check1=$ent;
    $check2=$aa;
    $prev_pdbchn = $pdbchn;

  } #foreach $line (@coord) 

  #----------------------------------------------------------------------------------
  # 3. Run Reduce_knot(): Determine if configuration is knotted, knot type and size |
  #----------------------------------------------------------------------------------

  foreach my $conf (@configs){


#     open(KNOT_FILE,">>$TMP/knot.dat");
#     print KNOT_FILE "$protein_array[$counter]\n"; 
#     close(KNOT_FILE);
      #`/bin/cp -f $TMP/k.$i $TMP/conf_in`;
      #print `$CGIPATH/Reduce_knot20`;
      print "$conf->{chn} $conf->{start} ";
      my @arg = $conf->{file};
      @arg = ("-r","$redknotbase.$conf->{chn}-$conf->{start}",@arg)
	  if defined $redknotbase;
      system "./Reduce_knot20", @arg;
      unlink $conf->{file};
  }

  #----------------------
  # 4. Delete old files |
  #----------------------

#    `rm -rf $TMP/k.* $TMP/conf_in`;

} #end OUTER LOOP -> proceed with next mmCIF
print "END\n";
