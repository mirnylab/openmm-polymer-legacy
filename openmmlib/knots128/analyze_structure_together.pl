#! /usr/bin/perl

`rm knot.dat`;

#----------------------------------------------------------------------------------------
# OUTER LOOP: Aquire list of proteins to be analyzed and loop over all mmCIF structures |
#----------------------------------------------------------------------------------------

@protein_array=`cat complete_list_2006.txt`;

$ending_value = scalar(@protein_array);
#$ending_value = 1;
for($counter=0 ; $counter < $ending_value ; $counter++)
{

  #----------------------------------------
  # 1. Get C-alpha coordinates from mmCIF |
  #----------------------------------------

  @coord = `awk '\$1==\"ATOM\" && \$4==\"CA\" {print \$5,\$8,\$9,\$11,\$12,\$13} END {print 0,0,0,0,0,0}' $protein_array[$counter]`;

  # $het := "$5"                 information regarding microheterogeneity (.,A,B,...)
  # $ent := "$8"                 signals when a new structure starts (1,2,...)
  # $aa  := "$9"                 sequence number of amino acid
  # $x   := "$11","$12","$13"    x,y,z coordinates

  #------------------------------------------------  
  # 2. Create configuration files for Reduce_knot |
  #------------------------------------------------

  $monomer_ctr       = 0; #atom number in structure 
  $configuration_ctr = 0; #configuration file number (k.0,k.1,...)
  $start             = 0;

  foreach $line (@coord) {
    #go through line by line
    ($het[$monomer_ctr],$ent,$aa,$x[$monomer_ctr],$y[$monomer_ctr],$z[$monomer_ctr]) = split(" ",$line,9999);

    #----------------------------------------------------------------------------------------------
    # 2.1. Split configuration if a new molecule starts or sequence of amino acids is interrupted |
    #      -> print to configuration file                                                         |
    #----------------------------------------------------------------------------------------------

    #if new section -> print out previous; else continue with next line
    #check1 and check2 store information of $ent and $aa from last iteration

    if( $monomer_ctr != 0           &&                          #not first line? AND
        ( $ent       != $check1 )   ||                          #( new molecule? OR
           ($aa < $check2)) { 
#          ($aa       != $check2+1   &&   $aa != $check2) ) ) {  #  interrupted sequence of amino acids) 


      #---------------------------------------------------------------------------------------------
      # 2.1.1. Check for microheterogeneity (file contains multiple coordinates for the same atom) |
      #---------------------------------------------------------------------------------------------

      undef $micro_het_symbol_number;
      undef @micro_het_symbol;
      $micro_het_symbol_number = 1;
      $micro_het_symbol[$micro_het_symbol_number] = '.'; 

      for ($i=$start,$j=1;$i<$monomer_ctr;$i++){
        #symbol already in symbol array ?
        for ($j=1,$k=0;$j<=$micro_het_symbol_number;$j++) {
           if($het[$i] eq $micro_het_symbol[$j]) {$k=1;}          #YES: set k=1
        }
        #no (k!=1): add symbol to symbol array; example: $micro_het_symbol[2]=A 
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
        $number_of_monomers=0;
        for ($i=$start,$j=1;$i<$monomer_ctr;$i++){
          $number_of_monomers++;
        }

        open(CONFIG_FILE,">k.$configuration_ctr");
        print CONFIG_FILE "t=0\n\n$number_of_monomers\n";

        #print coordinates
        for ($i=$start,$j=1;$i<$monomer_ctr;$i++){
          print CONFIG_FILE "$j\t$x[$i]\t$y[$i]\t$z[$i]\n";
          $j++;
        }
   
        close(CONFIG_FILE);
        $configuration_ctr++;
      } #end if($micro ...

      #----------------------------------------------------------------------------------------------
      # 2.1.2.b. Microheterogeneity: Print coordinates - each case gets separate configuration file |
      #----------------------------------------------------------------------------------------------

      else{
        for($k=2;$k<=$micro_het_symbol_number;$k++) {

          #determine number of monomers in configuration file
          $number_of_monomers=0; 
          for ($i=$start,$j=1;$i<$monomer_ctr;$i++)  {
            if ($het[$i] eq '.' || $het[$i] eq $micro_het_symbol[$k]) {
              $number_of_monomers++;
            }
          }

          open(CONFIG_FILE,">k.$configuration_ctr");
          print CONFIG_FILE "t=0\n\n$number_of_monomers\n";

          #print coordinates
          for ($i=$start,$j=1;$i<$monomer_ctr;$i++)  {
            if ($het[$i] eq '.' || $het[$i] eq $micro_het_symbol[$k]) {
              print CONFIG_FILE "$j\t$x[$i]\t$y[$i]\t$z[$i]\n"; 
              $j++;
            }
          }
        close(CONFIG_FILE);
        $configuration_ctr++;
        }
      } #end else{ ...

      $start=$monomer_ctr;

    } #end if( $monomer_ctr != 0 ... -> continue with next protein or protein part in mmCIF

    $monomer_ctr++;
    $check1=$ent;
    $check2=$aa;

  } #foreach $line (@coord) 

  #----------------------------------------------------------------------------------
  # 3. Run Reduce_knot(): Determine if configuration is knotted, knot type and size |
  #----------------------------------------------------------------------------------

  for($i=0;$i<$configuration_ctr;$i++) {

    open(KNOT_FILE,">>knot.dat");
    print KNOT_FILE "$protein_array[$counter]"; 
    close(KNOT_FILE);
    `cp k.$i conf_in`;
    `./Reduce_knot20 >> knot.dat`;
  }

  #----------------------
  # 4. Delete old files |
  #----------------------

  `rm -rf k.* conf_in`;

} #end OUTER LOOP -> proceed with next mmCIF
