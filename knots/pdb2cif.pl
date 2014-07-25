#!/usr/local/bin/perl
eval 'exec perl -S $0 "$@"'
    if $running_under_some_shell;
			# this emulates #! processing on NIH machines.
			# (remove #! line above if indigestible)

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;
			# process any FOO=bar switches

#
#
#
#
######################################################################
#
#
#                               pdb2cif.pl
#
#                          produced from pdb2cif.m4
#                         version 2.4.2   7 Oct 2004
#
#    a m4 macro program which produces pdb2cif.pl, pdb2cif.awk, pdb2cif.oawk
#
#            Scripts to filter a PDB entry and produce a CIF file.
#
#                       Phil Bourne (bourne@sdsc.edu)
#
#                  adapted to 6 Oct 95 cifdic.m95 0.7.28
#                  and later to 1 Jan 97 cif_mm.dic 0.9.01
#
#                                   by
#                           Herbert J. Bernstein
#               Bernstein+Sons, P.O. Box 177, Bellport, NY 11713
#           phone: 1-631-286-1339, email: yaya@bernstein-plus-sons.com
#                                  and
#                           Frances C. Bernstein
#               Bernstein+Sons, P.O. Box 177, Bellport, NY 11713
#           phone: 1-631-286-1339, email: fcb@bernstein-plus-sons.com
#
#  This work was supported in part by IUCr (for HJB), US NSF, PHS, NIH, 
#  NCRR, NIGMS, NLM and DOE (for FCB prior to 1998), US NSFgrant no. 
#  See H. Bernstein, F. Bernstein, P. E. Bourne "CIF Applications. VIII. 
#  pdb2cif: Translating PDB Entries into mmCIF Format", J. Appl. Cryst.,
#  31, pp. 282-295, 1998.
#
#**************************************************************************
#         THE CONVERSION FROM PDB FORMAT TO CIF FORMAT IS COMPLEX
#                     ******* USE WITH CAUTION *******
#                   COMMENTS AND SUGGESTIONS APPRECIATED
#      If you like the basic approach, thank Phil Bourne.  He did
#      the real work of creating pdb2cif.  If you have problems with
#      the adaptation to cif_mm.dic, tell yaya@bernstein-plus-sons.com
#
#
###########################################################################
#**************************************************************************
#
#
#  This version available via http from:
#
#      http://www.bernstein-plus-sons.com/software/pdb2cif
#      http://www.iucr.org/iucr-top/cif/software/pdb2cif
#            and the mirror sites of the IUCr
#      http://www.sdsc.edu/pb/pdb2cif/pdb2cif
#      http://ndbserver.rutgers.edu/NDB/mmcif/software/pdb2cif
#            and the mirror sites of the NDB
#
#  See the file README for instructions on use and installation
#
###########################################################################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
##############################################################################
#
# Version History:
#         See the file  CHANGES
#
#######################################################################

$[ = 1;			# set array base to 1

$comma = ',';
$lcaz = 'abcdefghijklmnopqrstuvwxyz';
$UCAZ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
$version = '2.4.2';
$version_date = ' 7 Oct 2004';
printf (("\n"));
printf (("###################################################\n"));
printf (("#                                                 #\n"));
printf (("#   Converted from PDB format to CIF format by    #\n"));
printf "#   pdb2cif version %-15s %11s   #\n", $version, $version_date;
printf (("#                       by                        #\n"));
printf (("# P.E. Bourne, H.J. Bernstein and F.C. Bernstein  #\n"));
printf (("#                                                 #\n"));
printf (("# http://www.bernstein-plus-sons.com/software     #\n"));
printf (("#                     /pdb2cif                    #\n"));
printf (("#   *** See the remarks at the end of this  ***   #\n"));
printf (("#   *** file for information on conversion  ***   #\n"));
printf (("#   *** of this entry and on the program    ***   #\n"));
printf (("#   ***               pdb2cif               ***   #\n"));
$my_at = '@';
printf (("# Please report problems to:                      #\n"));
printf ((('#         pdb2cif' . $my_at .

  "bernstein-plus-sons.com         #\n")));
printf (("###################################################\n\n\n"));
#
#  Set starting variables
#
#       The following flag is used to produce a more complete CIF entry,
#       i.e. data items are given, but with the value "?".
#       If you desire only the minimum set of data items comment out the
#       following one line:

#       verbose = "yes"
#
#       The following flag controls conversion of text fields using
#       the type-setting codes used in some PDB entries

#       convtext = "yes"
#
#       The following flag controls conversion of author and editor
#       names, "yes" to always convert according to the 1992 format
#       description, "conditional" to be controlled by convtext
$auth_convtext = 'yes';
#
#       uncomment the next line if convtext control of typesetting desired
#       auth_convtext = "conditional"

#
#       The following flag controls the distribution of entity_seq_num
#       to all atom site lines, uncomment if you do _not_ want
#       this distribution done, but want denser atom lists

#       dense_list = "yes"
#

#
#       The following flag controls the handling of Junior, Senior, etc
#       As of this writing, the mmCIF convention is to keep dynastic
#       indicators with the last name
#
$junior_on_last = 'yes';
#
#       The following flag controls the printing of TER records
#       The possible values are "yes" to print them, "no" to
#       suppress them, or "comment" to print them as comments
#
$print_ter = 'comment';
#
#       The following flag controls the use of the pdb2cif prefix
#       for tags which are not part of the mmCIF dictionary.
#       The possible values are "yes" to include the prefix or
#       "no" to suppress them
#
$use_pdb2cif_prefix = 'yes';
#
#
#       each of these variables may be reset from within the text
#       with a line of the form
#
#       #define variable value
#
#       e.g.   #define verbose yes
#       or     #define convtext no
#
#       In addition, to allow control of translation of author names
#       the variable "name" may be defined multiple times
#       with two values, in the form
#
#       #define name PDB_form name_value
#
#       where the PDB_form is the form of the name expected in the PDB
#       and name_value is the form to be used by this program.  All blanks
#       in either form must be replaced by "_".  For example, you can
#       give the following
#
#       #define name E.F.MEYER_JUNIOR Meyer Junior,_E.F.
#
#       If the same name is defined multiple times, only the last
#       translation given will be used.  The PDB_form is not case-sensitive,
#       but the name_value is.
#
#       Normally, the compliance_level for a PDB dataset should be
#       obtained from REMARK 4.  However, to facilitate processing
#       of pseudo-pdb datasets from non-PDB sources, the compliance
#       level may be set by
#
#       #define compliance_level level
#
#       where 2.0 is the only meaningful case of level at this time
#
#       If REMARK 4 contains a compliance level, that will apply in
#       parsing from that point onwards

$aniso_flag = 0;
$atom_align_ids = 'yes';
$atom_unique_ids = 'yes';
$atom_alt_flag = 1;
$atom_flag = 1;
$atom_flag_1 = 1;
$atom_flag_2 = 1;
$atom_ids_in_models = 0;
$atom_res_flag = 0;
$audit_flag = 0;
$author_flag = 1;
$bcid = '*';# character to use for blank chain id
$cispep_flag = 0;
$cit_flag = 0;# 1 is the primary citation
$compnd_flag = 1;
$compliance_level = 0.0;
$connect_flag = 0;
$conect_flag = 1;
$conect_flag_2 = 0;
$conect_id = 1;
$dbref_flag = 0;# count of dbref records
$end_flag = 0;
$entity_flag = 0;
$entity_mon_flag = 1;
$ent_non_poly_point{' '} = '';
$ent_non_poly_num{' '} = 0;
$ent_poly_point{' '} = '';
$ent_poly_num{' '} = 0;
$entity_seq_num_flag = 0;
$flush_ref = 0;
for ($X = 1; $X <= 999; ++$X) {
    $ftnote_flag{$X} = 0;
}
$ftnote_flag_old = -1;
$foot_flag = 0;
$formul_flag = 1;
$head_PDB_code = '.';
$helix_flag = 1;
$het_flag = 1;
$hetnam_flag = 1;
$hetsyn_flag = 1;
$hydbnd_flag = 0;
$jrnl_flag = 1;
$keywrd_flag = 0;
$link_flag = 0;
$mon_flag = 1;
$model_count = 0;
$model{$model_count} = '.';
$model_asn_low{$model_count} = 0;
$model_asn_high{$model_count} = 0;
$model_an_low{$model_count} = 0;
$model_an_high{$model_count} = 0;
$model_align_ids = 'yes';
$model_flag = '.';
$model_flags = 'no';
$modres_flag = 0;
$mtrix_flag = 0;
$nmr_flag = 0;
$nonp_flag = 1;
$num_non_poly_ents = 0;
$num_poly_ents = 0;
$num_res_name = 0;
$num_res_pair = 0;
$origx_flag = 0;
$pdb2cif_prefix = 'pdb2cif_';
if ($use_pdb2cif_prefix eq 'no') {
    $pdb2cif_prefix = '';
}
$record_number = ' ';
$remark_flag = 0;
$remark_header_flag = 0;
$revdat_flag = 1;
$res_flag = 0;
$res_res_flag = 0;
$scale_flag = 0;
$seqadv_flag = 0;# count of seqadv records
$seqres_flag = 1;
$sheet_flag = 0;
$sheet_flag_2 = 1;
$sigatm_flag = 0;
$siguij_flag = 0;
$site_flag = 1;
$sltbrg_flag = 0;
$s_o_flag = 0;
$ss_flag = 1;# tracks HELIX, TURN and SHEET
$ss_flag_2 = 1;
$ssbond_flag = 1;
$ter_flag = 0;
$turn_flag = 1;
$turn_flag_2 = 1;
$tvect_flag = 0;
$vol_flag = 0;
$warning_flag = 0;# count of warnings in warning list
$xlat_flag = 0;
$xlat_save = 0;

$flag = 0;# flags to correctly set ;
$previous_keyword = ' ';# prior to new keyword
$remark_number = 0;
$remark_number_old = 0;
$all_remarks = 0;
#
#       set up connect types
#
{
    $num_ctypes = (@connect_types = split(' ',

      '. . . . hydrog hydrog saltbr hydrog hydrog saltbr', 9999));
    if ($num_ctypes >= 1 && $connect_types[$num_ctypes] eq '' && ' ' eq ' ') {
	--$num_ctypes;
    }
}
#
#       set up conversion strings for residues
#
$numl = '0123456789';
$charl = $lcaz;
$charu = $UCAZ;
$chars = "+_*/!#\$,.;:?|{}()";
$charx = ($charl . $charu . $numl . $chars);
#       Define date format conversion arrays
#
#       mmm2mm[month_name] = month_ordinal
#       yyyy[2_digit_year] = 4_digit_year
#
$mmm2mm{'JAN'} = '01';
$mmm2mm{'FEB'} = '02';
$mmm2mm{'MAR'} = '03';
$mmm2mm{'APR'} = '04';
$mmm2mm{'MAY'} = '05';
$mmm2mm{'JUN'} = '06';
$mmm2mm{'JUL'} = '07';
$mmm2mm{'AUG'} = '08';
$mmm2mm{'SEP'} = '09';
$mmm2mm{'OCT'} = '10';
$mmm2mm{'NOV'} = '11';
$mmm2mm{'DEC'} = '12';
$mmm2mm{'Jan'} = '01';
$mmm2mm{'Feb'} = '02';
$mmm2mm{'Mar'} = '03';
$mmm2mm{'Apr'} = '04';
$mmm2mm{'May'} = '05';
$mmm2mm{'Jun'} = '06';
$mmm2mm{'Jul'} = '07';
$mmm2mm{'Aug'} = '08';
$mmm2mm{'Sep'} = '09';
$mmm2mm{'Oct'} = '10';
$mmm2mm{'Nov'} = '11';
$mmm2mm{'Dec'} = '12';
#
for ($yy = 0; $yy < 100; ++$yy) {
    $yyyy{$yy + 0} = $yy + 1900;
    if ($yy < 70) {
	$yyyy{$yy + 0} += 100;
    }
}
#
#       Define lists of amino acids and nucleic acids
{
    $num_aa = (@aa_list = split(' ',

      ('ABU ACD ALA ALB ALI ARG ARO ASN ASP ASX' .

      ' BAS CYS GLN GLU GLX GLY HIS HYP ILE LEU' .

      ' LYS MET PCA PHE PRO SER THR TRP TYR VAL'), 9999));
    if ($num_aa >= 1 && $aa_list[$num_aa] eq '' && ' ' eq ' ') {
	--$num_aa;
    }
}
{
    $num_na = (@na_list = split(' ', ('A +A C +C G +G I +I T +T U +U'),

      9999));
    if ($num_na >= 1 && $na_list[$num_na] eq '' && ' ' eq ' ') {
	--$num_na;
    }
}

#
#  Formulae and Molecular Weights for Standard Residues from
#  the 1992 PDB format Description
#  AARESES has the Amino Acids and Miscellaneous Residues
#  NARESES has the Nucleotides
#
# Name; Code;Formula;Mol. Wt.
#

{
    $num_AARESES = (@AARESES_list = split(/\|/,

      ('Alanine;ALA;C3 H7 N1 O2;89.09|' . 'Arginine;ARG;C6 H14 N4 O2;174.20|'

      . 'Asparagine;ASN;C4 H8 N2 O3;132.12|' .

      'Aspartic acid;ASP;C4 H7 N1 O4;133.10|' .

      'ASP/ASN ambiguous;ASX;C4 H7.5 N1.5;132.61|' .

      'Cysteine;CYS;C3 H7 N1 O2 S1;121.15|' .

      'Glutamine;GLN;C5 H10 N2 O3;146.15|' .

      'Glutamic acid;GLU;C5 H9 N1 O4;147.13|' .

      'GLU/GLN ambiguous;GLX;C5 H9.5 N1.5 O3.5;146.64|' .

      'Glycine;GLY;C2 H5 N1 O2;75.07|' . 'Histidine;HIS;C6 H9 N3 O2;155.16|' .

      'Isoleucine;ILE;C6 H13 N1 O2;131.17|' .

      'Leucine;LEU;C6 H13 N1 O2;131.17|' . 'Lysine;LYS;C6 H14 N2 O2;146.19|' .

      'Methionine;MET;C5 H11 N1 O2 S1;149.21|' .

      'Phenylalanine;PHE;C9 H11 N1 O2;165.19|' .

      'Proline;PRO;C5 H9 N1 O2;115.13|' . 'Serine;SER;C3 H7 N1 O3;105.09|' .

      'Threonine;THR;C4 H9 N1 O3;119.12|' .

      'Tryptophan;TRP;C11 H12 N2 O2;204.23|' .

      'Tyrosine;TYR;C9 H11 N1 O3;181.19|' . 'Valine;VAL;C5 H11 N1 O2;117.15|'

      . 'Undetermined;UNK;C5 H6 N1 O3;128.16|' .

      'Acetic Acid;ACE;C2 H4 O2;60.05|' . 'Formic Acid;FOR;C1 H2 O2;40.03|' .

      'Water;HOH;H2 O1;18.015'), 9999));
    if ($num_AARESES >= 1 && $AARESES_list[$num_AARESES] eq '' &&

      '|' eq ' ') {
	--$num_AARESES;
    }
}

for ($naa = 1; $naa <= $num_AARESES; ++$naa) {
    {
	$nxx = (@naa_split = split(/;/, $AARESES_list[$naa], 9999));
	if ($nxx >= 1 && $naa_split[$nxx] eq '' && ';' eq ' ') {
	    --$nxx;
	}
    }
    $res_id{$naa} = $naa_split[2];
    $res_count{$naa_split[2]} = 0;
    $res_name{$naa_split[2]} = $naa_split[1];
    $res_formul{$naa_split[2]} = $naa_split[3];
}

{
    $num_NARESES = (@NARESES_list = split(/\|/,

      ('Adenosine;  A;C10 H14 N5 O7 P1;347.22|' .

      'Modified Adenosine; +A;.;347.22|' .

      '1-Methyladenosine;1MA;C11 H16 N5 O7 P1;361.25|' .

      'Cytidine;  C;C9 H14 N3 O8 P1;323.20|' .

      'Modified Cytidine; +C;.;323.20|' .

      '5-Methylcytidine;5MC;C10 H16 N3 O8 P1;337.23|' .

      "2'-O-Methylcytidine;OMC;C10 H17 N3 O8 P1;338.23|" .

      'Guanosine;  G;C10 H14 N5 O8 P1;363.22|' .

      'Modified Guanosine; +G;.;363.22|' .

      '1-Methylguanosine;1MG;C11 H16 N5 O8 P1;377.25|' .

      'N2-Methylguanosine;2MG;C11 H16 N5 O8 P1;377.25|' .

      'N2-Dimethylguanosine;M2G;C12 H18 N5 O8 P1;391.28|' .

      '7-Methylguanosine;7MG;C11 H10 N5 O8 P1;377.25|' .

      "2'-O-Methylguanosine;OMG;C11 H16 N5 O8 P1;377.25|" .

      'Wybutosine; YG;C21 H26 N6 O11 P1;587.48|' .

      'Inosine;  I;C10 H13 N4 O8 P1;348.21|' .

      'Modified Inosine; +I;.;348.21|' .

      'Thymidine;  T;C10 H15 N2 O8 P1;322.21|' .

      'Modified Thymidine; +T;.;322.21|' .

      'Uridine;  U;C9 H13 N2 O9 P1;324.18|' . 'Modified Uridine; +U;.;324.18|'

      . 'Dihydrouridine;H2U;C9 H15 N2 O9 P1;326.20|' .

      'Ribosylthymidine;5MU;C10 H16 N2 O10 P1;355.22|' .

      'Pseudouridine;PSU;C9 H13 N2 O9 P1;324.18|'), 9999));
    if ($num_NARESES >= 1 && $NARESES_list[$num_NARESES] eq '' &&

      '|' eq ' ') {
	--$num_NARESES;
    }
}
for ($nna = 1; $nna <= $num_NARESES; ++$nna) {
    {
	$nxx = (@nna_split = split(/;/, $NARESES_list[$nna], 9999));
	if ($nxx >= 1 && $nna_split[$nxx] eq '' && ';' eq ' ') {
	    --$nxx;
	}
    }
    $res_id{$num_AARESES + $nna} = $nna_split[2];
    $res_count{$nna_split[2]} = 0;
    $res_name{$nna_split[2]} = $nna_split[1];
    $res_formul{$nna_split[2]} = $nna_split[3];
}
$num_OTRESES = 0;

#
#       Element Lists to check atom types
#
{
    $num_per_tab = (@periodic_table = split(' ',

      (' . D  ' . 'H                                                  HE ' .

      'LI BE                               B  C  N  O  F  NE ' .

      'NA MG                               AL SI P  S  CL AR ' .

      'K  CA SC TI V  CR MN FE CO NI CU ZN GA GE AS SE BR KR ' .

      'RB SR Y  ZR NB MO TC RU RH PD AG CD IN SN SB TE I  XE ' . 'CS BA ' .

      '      LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU ' .

      '         HF TA W  RE OS IR PT AU HG TL PB BI PO AT RN ' . 'FR RA ' .

      '      AC TH PA U  NP PU AM CM BK CF ES FM MD NO LR ' .

      '         KU HA SG NS HS' . '         DB JL RF BH HN MT'), 9999));
    if ($num_per_tab >= 1 && $periodic_table[$num_per_tab] eq '' &&

      ' ' eq ' ') {
	--$num_per_tab;
    }
}
{
    $num_aa_na_el = (@standard_res_elements = split(' ', ('. C H D N O P S'),

      9999));
    if ($num_aa_na_el >= 1 && $standard_res_elements[$num_aa_na_el] eq '' &&

      ' ' eq ' ') {
	--$num_aa_na_el;
    }
}
{
    $num_one_let_el = (@one_letter_elements = split(' ',

      ('. B C H D N O F P S K V Y I W'), 9999));
    if ($num_one_let_el >= 1 && $one_letter_elements[$num_one_let_el] eq '' &&

      ' ' eq ' ') {
	--$num_one_let_el;
    }
}

#
#       Define special name suffixes to move away from family names
#
{
    $num_suffix = (@suffix_list = split(' ',

      ('JUNIOR SENIOR JR SR JR. SR.' .

      ' I II III IV V VI VII VIII IX X XI XII'), 9999));
    if ($num_suffix >= 1 && $suffix_list[$num_suffix] eq '' && ' ' eq ' ') {
	--$num_suffix;
    }
}
{
    $xnum_suffix = (@rep_suffix_list = split(' ',

      ('Junior Senior Junior Senior Junior Senior' .

      ' I II III IV V VI VII VIII IX X XI XII'), 9999));
    if ($xnum_suffix >= 1 && $rep_suffix_list[$xnum_suffix] eq '' &&

      ' ' eq ' ') {
	--$xnum_suffix;
    }
}
for ($i = 1; $i <= $xnum_suffix; ++$i) {
    $rep_suffix{$suffix_list[$i]} = $rep_suffix_list[$i];
}

#
#        Setup charge conversions
#

#        No Charge
$charge{'0 '} = '  ';
$charge{'00'} = '  ';
$charge{'  '} = '  ';

#        1+
$charge{'1+'} = '1+';
$charge{'+1'} = '1+';
$charge{'I '} = '1+';
$charge{'i '} = '1+';
$charge{'1 '} = '1+';
$charge{'+ '} = '1+';

#        1-
$charge{'1-'} = '1-';
$charge{'-1'} = '1-';
$charge{'- '} = '1-';

#        2+
$charge{'2+'} = '2+';
$charge{'+2'} = '2+';
$charge{'II'} = '2+';
$charge{'ii'} = '2+';
$charge{'2 '} = '2+';
$charge{'++'} = '2+';

#        2-
$charge{'2-'} = '2-';
$charge{'-2'} = '2-';
$charge{'--'} = '2-';

#        3+
$charge{'3+'} = '3+';
$charge{'+3'} = '3+';
$charge{'3 '} = '3+';

#        3-
$charge{'3-'} = '3-';
$charge{'-3'} = '3-';

#        4+
$charge{'4+'} = '4+';
$charge{'+4'} = '4+';
$charge{'4 '} = '4+';

#        4-
$charge{'4-'} = '4-';
$charge{'-4'} = '4-';

#        5+
$charge{'5+'} = '5+';
$charge{'+5'} = '5+';
$charge{'5 '} = '5+';

#        5-
$charge{'5-'} = '5-';
$charge{'-5'} = '5-';
$charge{'5 '} = '5-';

#        6+
$charge{'6+'} = '6+';
$charge{'+6'} = '6+';
$charge{'6 '} = '6+';

#        6-
$charge{'6-'} = '6-';
$charge{'-6'} = '6-';

#        7+
$charge{'7+'} = '7+';
$charge{'+7'} = '7+';
$charge{'7 '} = '7+';

#        7-
$charge{'7-'} = '7-';
$charge{'-7'} = '7-';

#        8+
$charge{'8+'} = '8+';
$charge{'+8'} = '8+';
$charge{'8 '} = '8+';

#        8-
$charge{'8-'} = '8-';
$charge{'-8'} = '8-';

#        9+
$charge{'9+'} = '9+';
$charge{'+9'} = '9+';
$charge{'9 '} = '9+';

#        9-
$charge{'9-'} = '9-';
$charge{'-9'} = '9-';

# End of BEGIN statement
#
#
#   Flag all lines as untranslated unless proven otherwise.
#
#   Process #define (or #def)
#
#   Determine whether this is a new keyword, if so and flag is set
#   terminate free text with a ;  Also discard noise lines less than
#   6 characters long, and pad other lines to 80 characters with blanks
#
#   Ensure that the record name used is separated from following info
#

while (<>) {
    chop;	# strip record separator
    @Fld = split(' ', $_, 9999);
    if ($end_flag == 0) {
	$xlat_save = $xlat_flag;
	$non_xlated{++$xlat_flag} = $_;
	$first_field = $Fld[1];
	if ($#Fld > 1 && ($Fld[1] eq '#def' || $Fld[1] eq '#define')) {
	    {
		$lx_tl = length($Fld[2]);
		$tx_tl = $Fld[2];
		$var_name = '';
		for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
		    $cx_tl = substr($tx_tl, $ix_tl, 1);
		    $cx_tl = substr(($lcaz . $cx_tl), index(($UCAZ . $cx_tl),

		      $cx_tl), 1);
		    $var_name = ($var_name . $cx_tl);
		}
	    }

	    $var_value = '';
	    if ($#Fld > 2) {
		$lx_tl = length($Fld[3]);
		$tx_tl = $Fld[3];
		$var_value = '';
		for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
		    $cx_tl = substr($tx_tl, $ix_tl, 1);
		    $cx_tl = substr(($lcaz . $cx_tl), index(($UCAZ . $cx_tl),

		      $cx_tl), 1);
		    $var_value = ($var_value . $cx_tl);
		}
	    }

	    if ($var_name eq 'verbose' && ($var_value eq 'yes' ||

	      $var_value eq 'no')) {
		$verbose = $var_value;
		$xlat_flag = $xlat_save;
	    }
	    else {
		if ($var_name eq 'convtext' && ($var_value eq 'yes' ||

		  $var_value eq 'no')) {
		    $convtext = $var_value;
		    $xlat_flag = $xlat_save;
		}
		else {
		    if ($var_name eq 'auth_convtext' && 
		    ($var_value eq 'yes' || $var_value eq 'no' ||

		      $var_value eq 'conditional')) {
			$auth_convtext = $var_value;
			$xlat_flag = $xlat_save;
		    }
		    else {
			if ($var_name eq 'dense_list' && 
			($var_value eq 'yes' || $var_value eq 'no')) {
			    $dense_list = $var_value;
			    $xlat_flag = $xlat_save;
			}
			else {
			    if ($var_name eq 'junior_on_last' && 
			    ($var_value eq 'yes' || $var_value eq 'no')) {
				$junior_on_last = $var_value;
				$xlat_flag = $xlat_save;
			    }
			    else {
				if ($var_name eq 'print_ter' && 
				($var_value eq 'yes' || $var_value eq 'no' ||

				  $var_value eq 'comment')) {
				    $print_ter = $var_value;
				    $xlat_flag = $xlat_save;
				}
				else {
				    if ($var_name eq 'compliance_level') {
					$compliance_level = $var_value;
					$xlat_flag = $xlat_save;
				    }
				    else {
					if ($var_name eq 'use_pdb2cif_prefix')

					  {
					    if ($var_value eq 'no') {
						;
					    }
					    else {
						if ($var_value eq 'yes') {
						    $use_pdb2cif_prefix =

						      'yes';
						    $pdb2cif_prefix =

						      'pdb2cif_';
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	    if ($#Fld > 3 && $var_name eq 'name') {
		$ll = length($Fld[3]);
		$PDB_form = '';
		for ($i = 1; $i <= $ll; ++$i) {
		    $cc = substr($Fld[3], $i, 1);
		    if ($cc eq '_') {	#???
			$cc = ' ';
		    }
		    $PDB_form = ($PDB_form . $cc);
		}
		$ll = length($Fld[4]);
		$name_value = '';
		for ($i = 1; $i <= $ll; ++$i) {
		    $cc = substr($Fld[4], $i, 1);
		    if ($cc eq '_') {	#???
			$cc = ' ';
		    }
		    $name_value = ($name_value . $cc);
		}
		{
		    #
		    # apply PDB typsetting codes if any to a line
		    #
		    {
			$lx_tl = length($PDB_form);
			$tx_tl = $PDB_form;
			$lostr = '';
			for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
			    $cx_tl = substr($tx_tl, $ix_tl, 1);
			    $cx_tl = substr(($lcaz . $cx_tl),

			      index(($UCAZ . $cx_tl), $cx_tl), 1);
			    $lostr = ($lostr . $cx_tl);
			}
		    }

		    $lstr = length($lostr);
		    $mystr = '';
		    $pchar = ' ';
		    for ($qtsi = 1; $qtsi <= $lstr; ++$qtsi) {
			$mychar = substr($lostr, $qtsi, 1);
			if ($pchar eq ' ' || $pchar eq ',' || $pchar eq '.' ||

			  $pchar eq '(' || $pchar eq '*' || $pchar eq '/') {
			    $mychar = substr(($UCAZ . $mychar),

			      index(($lcaz . $mychar), $mychar), 1);
			}
			if (($mychar ne '*' && $mychar ne "\$" &&

			  $mychar ne '/') || ($mychar eq $pchar)) {	#???
			    $mystr = ($mystr . $mychar);
			}
			if ($pchar eq '/') {
			    if ($mychar eq "\$" || $mychar eq '-') {
				$pchar = $mychar;
			    }
			}
			else {
			    $pchar = $mychar;
			}
		    }
		    $ret_val = $mystr;

		    $PDB_form = $ret_val;
		}

		{
		    $lx_tu = length($PDB_form);
		    $tx_tu = $PDB_form;
		    $PDB_form = '';
		    for ($ix_tu = 1; $ix_tu <= $lx_tu; ++$ix_tu) {
			$cx_tu = substr($tx_tu, $ix_tu, 1);
			$cx_tu = substr(($UCAZ . $cx_tu),

			  index(($lcaz . $cx_tu), $cx_tu), 1);
			$PDB_form = ($PDB_form . $cx_tu);
		    }
		}

		$rep_name{$PDB_form} = $name_value;
		$xlat_flag = $xlat_save;
	    }
	}
	if (length($_) > 5) {
	    if (length($_) < 80) {
		$_ = (($_) .

		  substr(('                                        ' .

		  '                                        '), 1,

		  80 - length($_)));
	    }
	    if (length($first_field) > 6) {
		$first_field = substr($Fld[1], 1, 6);
	    }
	    if ($first_field ne $previous_keyword && $flag ne '0') {	#???	#???
		printf (("; \n\n"));
		$flag = '0';
		$previous_keyword = $first_field;
	    }
	    else {
		$previous_keyword = $first_field;
	    }
	}
    }
    else {
	$_ = '';
	$first_field = '';
	if ($previous_keyword ne '' && $flag ne '0') {	#???
	    printf (("; \n\n"));
	    $flag = '0';
	}
	$previous_keyword = '';
    }

    #
    #  Print out any accumulated COMPND, SOURCE, TITLE or CAVEAT information

    if ($compnd_flag ne '1' && $first_field ne 'COMPND' &&	#???

      $first_field ne 'TITLE' && $first_field ne 'CAVEAT' &&

      $first_field ne 'SOURCE') {
	printf (("\n\n"));
	printf (("##################\n"));
	printf (("#                #\n"));
	printf (("#  STRUCT        #\n"));
	printf (("#                #\n"));
	printf (("##################\n\n"));
	printf (("loop_\n_struct.entry_id\n_struct.title\n"));
	printf "  %s\n", $head_PDB_code;
	printf "; %s\n", $compnd{1};
	for ($i = 2; $i < $compnd_flag; ++$i) {
	    printf "  %s\n", $compnd{$i};
	}
	printf (("; \n"));
	$compnd_flag = 1;
    }

    #=========================================================================
    #   Keyword ATOM or HETATM or TER
    #
    #  atom pdb type      [ 1- 6]  = _atom_site.group_PDB
    #  atom serial number [ 7-11]  = _atom_site.id
    #  atom type          [13-14]  = _atom_site.type_symbol
    #        (first 2 characters of atom name)
    #  atom name          [13-16]  = _atom_site.label_atom_id
    #  alternate location [17]     = _atom_site.label_alt_id
    #  residue name       [18-20]  = _atom_site.label_comp_id
    #  chain identifier   [22]     = _atom_site.label_asym_id
    #  residue seq no.    [23-26]  = _atom_site.auth_seq_id
    #  insertion code     [27]     =  appended to residue sequence no.
    #  x-coordinate       [31-38]  = _atom_site.cartn_x
    #  y-coordinate       [39-46]  = _atom_site.cartn_y
    #  z-coordinate       [47-54]  = _atom_site.cartn_z
    #  occupancy          [55-60]  = _atom_site.occupancy
    #  temperature factor [61-66]  = _atom_site.B_iso_or_equiv
    #  footnote number    [68-70]  = _atom_site.footnote_id
    #    (February 1992 PDB format)
    #  segment identifier [73-76]  = _atom_site.auth_asym_id
    #    (February 1996 PDB format)
    #  element symbol     [77-78]  = _atom_site.type_symbol
    #    (February 1996 PDB format)
    #  charge on atom     [79-80]  =   append to _atom_site.type_symbol
    #    (February 1996 PDB format)
    #
    # Information on non_standard monomers and non-polymers derived from
    # HET and FORMUL records is presented here using additional information
    # derived from ATOM and HETATM records.
    # The assignment of non-standard monomers versus non-polymers
    # is tricky and unlikely to be correct for all entries. Assignment is
    # based on the following rules:
    #   i) If the HET has a chain id then it must be non-standard (this
    #      is not complete since single chains do not have an chain id
    #      assigned.
    #   ii)If FORMUL places assigns a HET to a component number among
    #      the SEQRES components, the HET must be non-standard
    #
    #

    if ($first_field eq 'ATOM' || $first_field eq 'HETATM' ||

      $first_field eq 'TER') {
	$xlat_flag = $xlat_save;
	# parse field and save ATOM/HETATM/TER info
	# Since atoms are not necessarily numbered consecutively maintain
	# a complete conesecutive list 1 -> atom_flag and a partial
	# list for use by CONECT which references the atom_number
	#
	$atom_pdb = substr(($_), 7, 5);
	$atom_number{$atom_flag} = substr(($_), 7, 5);
	if ($atom_point{0 + $atom_number{$atom_flag}} eq '') {
	    $atom_point{0 + $atom_number{$atom_flag}} = $atom_flag;
	}
	else {
	    $atom_unique_ids = 'no';
	}
	$atom_name{$atom_flag} = substr(($_), 13, 4);
	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = $atom_name{$atom_flag};
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $temp_name = $ret_val;
	}
	$atom_name{$atom_flag} = $temp_name;
	$residue_name{$atom_flag} = substr(($_), 18, 3);
	$temp_name = substr($temp_name, 1, 2);
	{
	    $lx_tu = length($temp_name);
	    $tx_tu = $temp_name;
	    $temp_name = '';
	    for ($ix_tu = 1; $ix_tu <= $lx_tu; ++$ix_tu) {
		$cx_tu = substr($tx_tu, $ix_tu, 1);
		$cx_tu = substr(($UCAZ . $cx_tu), index(($lcaz . $cx_tu),

		  $cx_tu), 1);
		$temp_name = ($temp_name . $cx_tu);
	    }
	}

	if (index($UCAZ, substr($temp_name, 1, 1)) == 0) {
	    $temp_name = (' ' . substr($temp_name, 2, 1));
	}
	$xtemp_name = substr($temp_name, 1, 2);
	$ytemp_name = substr($temp_name, 2, 1);
	if (substr($xtemp_name, 1, 1) eq ' ') {
	    $xtemp_name = $ytemp_name;
	}
	$found = 'false';
	if ($first_field eq 'ATOM') {
	    ++$atom_res_flag;
	    if ($res_count{$residue_name{$atom_flag}} eq '') {
		$res_count{$residue_name{$atom_flag}} = 0;
		$res_formul{$residue_name{$atom_flag}} = '.';
		$res_name{$residue_name{$atom_flag}} = '.';
		++$num_OTRESES;
		$res_id{$num_AARESES + $num_NARESES + $num_OTRESES} =

		  $residue_name{$atom_flag};
		$warning_list{++$warning_flag} =

		  ('#=# ATOM_SITE: Residue name ' . $residue_name{$atom_flag}

		  . " not in standard residue list \n");
	    }
	    ++$res_count{$residue_name{$atom_flag}};
	    for ($ii = 1; $ii <= $num_aa_na_el && $found eq 'false'; ++$ii) {
		if ($xtemp_name eq $standard_res_elements[$ii]) {	#???
		    $found = 'true';
		}
	    }
	}
	else {
	    if (index($xhet_formula{$residue_name{$atom_flag}},

	      $xtemp_name) > 0) {
		$found = 'true';
	    }
	    else {
		if (index($xhet_formula{$residue_name{$atom_flag}},

		  (' ' . $ytemp_name . ' ')) > 0) {
		    $found = 'true';
		    if ($het_conv{($residue_name{$atom_flag} . '|' .

		      $xtemp_name . '|' . $ytemp_name)} eq '') {
			$warning_list{++$warning_flag} =

			  ('#=# ATOM_SITE: Het group ' .

			  $residue_name{$atom_flag} . '; atom type ' .

			  $xtemp_name . ' converted to ' . $ytemp_name .

			  "\n");
			$het_conv{($residue_name{$atom_flag} . '|' .

			  $xtemp_name . '|' . $ytemp_name)} = 'done';
		    }
		    $xtemp_name = $ytemp_name;
		    $temp_name = $ytemp_name;
		}
		else {
		    for ($ii = 1; $ii <= $num_per_tab && $found eq 'false';

		    ++$ii) {
			if ($xtemp_name eq $periodic_table[$ii]) {	#???
			    $found = 'true';
			}
		    }
		}
	    }
	}
	if ($found eq 'false') {
	    $temp_name = ' .';
	    if ($first_field eq 'ATOM') {
		for ($ii = 1; $ii <= $num_aa_na_el && $found eq 'false';

		++$ii) {
		    if ($ytemp_name eq $standard_res_elements[$ii]) {	#???
			$found = 'true';
		    }
		}
	    }
	    else {
		for ($ii = 1; $ii <= $num_one_let_el && $found eq 'false';

		++$ii) {
		    if ($ytemp_name eq $one_letter_elements[$ii]) {	#???
			$found = 'true';
		    }
		}
	    }
	    if ($found eq 'true') {
		$temp_name = (' ' . $ytemp_name);
	    }
	    $warning_list{++$warning_flag} = ('#=# ATOM_SITE: Site ' .

	      $atom_pdb . '; unexpected atom type ' . $xtemp_name .

	      ' converted to ' . $temp_name . "\n");
	}
	$atom_type{$atom_flag} = substr($temp_name, 1, 2);
	$atom_alt_location{$atom_flag} = substr(($_), 17, 1);
	$chain_id{$atom_flag} = substr(($_), 22, 1);
	$residue_seq_number{$atom_flag} = substr(($_), 23, 5);
	$atom_x{$atom_flag} = substr(($_), 31, 8);
	$atom_y{$atom_flag} = substr(($_), 39, 8);
	$atom_z{$atom_flag} = substr(($_), 47, 8);
	$atom_occ{$atom_flag} = substr(($_), 55, 6);
	$B_or_U{$atom_flag} = substr(($_), 61, 6);
	$footnote_number{$atom_flag} = substr(($_), 68, 3);
	if ($compliance_level >= 2.0) {
	    $atom_seg_id{$atom_flag} = substr(($_), 73, 4);
	    $atom_type{$atom_flag} = substr(($_), 77, 4);
	}
	{
	    #
	    #
	    # fix up atom_type (atom symbol and charge)
	    #
	    $temp_a_type = ($atom_type{$atom_flag} . '    ');
	    $orig_charge = substr($temp_a_type, 3, 2);
	    if ($orig_charge ne '  ') {
		if (substr($temp_a_type, 3, 1) eq ' ') {
		    $temp_a_type = (substr($temp_a_type, 1,

		      2) . substr($temp_a_type, 4, 1) . ' ');
		    $orig_charge = substr($temp_a_type, 3, 2);
		}
		$temp_charge = $charge{$orig_charge};
		if ($temp_charge ne '') {
		    $temp_a_type = (substr($temp_a_type, 1,

		      2) . $temp_charge);
		}
	    }
	    if (substr($temp_a_type, 2, 1) eq ' ') {
		$temp_a_type = (' ' . substr($temp_a_type, 1,

		  1) . substr($temp_a_type, 3, 2));
	    }
	    if ($temp_a_type eq '    ') {
		$temp_a_type = ' .  ';
	    }
	    if (substr(($temp_a_type . '  '), 3, 2) eq '  ') {
		$temp_a_type = substr($temp_a_type, 1, 2);
	    }
	    $ret_val = $temp_a_type;

	    $temp_type = $ret_val;
	}
	$atom_type{$atom_flag} = $temp_type;
	if ($atom_x{$atom_flag} eq '        ') {
	    $atom_x{$atom_flag} = '    .   ';
	}
	if ($atom_y{$atom_flag} eq '        ') {
	    $atom_y{$atom_flag} = '    .   ';
	}
	if ($atom_z{$atom_flag} eq '        ') {
	    $atom_z{$atom_flag} = '    .   ';
	}
	if ($atom_occ{$atom_flag} eq '      ') {
	    $atom_occ{$atom_flag} = '   .  ';
	}
	if ($B_or_U{$atom_flag} eq '      ') {
	    $B_or_U{$atom_flag} = '   .  ';
	}
	if ($atom_flag - $atom_number{$atom_flag} != 0) {
	    $atom_align_ids = 'no';
	}
	if ($model_flag ne '.') {
	    $model_asn_high{$model_count} = $atom_flag;
	    $model_an_high{$model_count} = $atom_number{$atom_flag};
	    if ($model_asn_low{$model_count} eq $model_asn_high{$model_count})	#???

	      {
		$model_an_low{$model_count} = $atom_number{$atom_flag};
	    }
	    if (($atom_flag - $model_asn_low{$model_count} + 1 -

	      $atom_number{$atom_flag}) != 0) {
		$model_align_ids = 'no';
	    }
	}
	#
	#
	#   flag atom as ATOM or HETATM or TER
	#
	if ($first_field eq 'ATOM') {
	    $atom_or_het{$atom_flag} = 'ATOM';
	}
	if ($first_field eq 'HETATM') {
	    $atom_or_het{$atom_flag} = 'HETATM';
	}
	if ($first_field eq 'TER') {
	    $atom_or_het{$atom_flag} = 'TER';
	    #
	    #   set alternate location value if blank
	    #
	    ;
	}
	if ($atom_alt_location{$atom_flag} eq ' ') {
	    $atom_alt_location{$atom_flag} = '.';

	    #
	    #   make a list of alternative atoms
	    #
	    ;
	}
	if ($atom_alt_location{$atom_flag} ne '.') {
	    $at_alt = $atom_alt_location{$atom_flag};
	    $atom_alt_list{$at_alt}++;
	}
	#
	#   set footnote value if blank
	#
	if ($footnote_number{$atom_flag} eq '   ') {
	    $footnote_number{$atom_flag} = ' . ';
	    #
	    #   set chain_id and entity_id to bcid for ATOM records if blank
	    #
	    ;
	}
	if (($first_field eq 'ATOM' ||

	  $first_field eq 'TER') && $chain_id{$atom_flag} eq ' ') {
	    $chain_id{$atom_flag} = $bcid;
	    $entity_id{$atom_flag} = $bcid;
	}
	#
	#   set chain_id to . for HETATM records if blank
	#
	if ($first_field eq 'HETATM' && $chain_id{$atom_flag} eq ' ') {
	    $chain_id{$atom_flag} = '.';
	    if ($num_poly_ents == 1 && $entities{1} eq $bcid) {	#???
		$chain_id{$atom_flag} = $bcid;
	    }
	}

	#
	#   set entity_id to chain_id for ATOM and TER records
	#
	if ($first_field eq 'ATOM' && $chain_id{$atom_flag} ne ' ') {
	    $entity_id{$atom_flag} = $chain_id{$atom_flag};
	}
	if ($first_field eq 'TER' && $chain_id{$atom_flag} ne ' ') {
	    $entity_id{$atom_flag} = $chain_id{$atom_flag};
	}

	#
	#  set _entity.id to residue_name for HETATM records
	#
	if ($first_field eq 'HETATM') {
	    $entity_id{$atom_flag} = $residue_name{$atom_flag};
	    $hetatm_entity = $residue_name{$atom_flag};
	    $ent_non_poly_id{$hetatm_entity}++;
	    if ($ent_non_poly_id{$hetatm_entity} == 1) {
		$next_non_poly_id = $ent_non_poly_point{' '};
		$prev_non_poly_id = ' ';
		while ($next_non_poly_id ne '') {
		    $prev_non_poly_id = $next_non_poly_id;
		    $next_non_poly_id =

		      $ent_non_poly_point{$prev_non_poly_id};
		}
		$ent_non_poly_point{$prev_non_poly_id} = $hetatm_entity;
		$ent_non_poly_point{$hetatm_entity} = '';
		++$num_non_poly_ents;
		$ent_non_poly_num{$hetatm_entity} = $num_non_poly_ents;
	    }
	    if ($entity_seq_num{$residue_name{$atom_flag}} ne '' &&

	      $entity_seq_num{$residue_name{$atom_flag}} + 0 <=

	      $num_poly_ents) {
		$entity_id{$atom_flag} = $chain_id{$atom_flag};
	    }
	}
	#
	#   define  _entities for polypeptide chains or DNA strands
	#   ie these are _entity_poly. Done by checking for chain in chain_id
	#   in ATOM records

	if ($first_field eq 'ATOM') {
	    $atom_entity = $chain_id{$atom_flag};
	    $ent_poly_id{$atom_entity}++;
	    if ($ent_poly_id{$atom_entity} == 1) {
		$next_poly_id = $ent_poly_point{' '};
		$prev_poly_id = ' ';
		while ($next_poly_id ne '') {
		    $prev_poly_id = $next_poly_id;
		    $next_poly_id = $ent_poly_point{$prev_poly_id};
		}
		$ent_poly_point{$prev_poly_id} = $atom_entity;
		$ent_poly_point{$atom_entity} = '';
		++$num_poly_ents;
		$ent_poly_num{$atom_entity} = $num_poly_ents;
		$entity_seq_num{$atom_entity} = $num_poly_ents;
		$entities{$num_poly_ents} = $atom_entity;
	    }
	}
	++$atom_flag;
    }

    #=====================================================================
    #   Keyword ANISOU
    #
    #
    #       atom serial number   =  matched via pointers to ATOM/HETATM
    #       atom type            =  dropped, taken from ATOM/HETATM
    #       atom name            =  dropped, taken from ATOM/HETATM
    #       alternate location   =  dropped, taken from ATOM/HETATM
    #       residue name         =  dropped, taken from ATOM/HETATM
    #       chain identifier     =  dropped, taken from ATOM/HETATM
    #       residue sequence no. =  dropped, taken from ATOM/HETATM
    #       insertion code       =  dropped, taken from ATOM/HETATM
    #
    #       
    #
    #     Note the different order
    #     PDB     CIF
    # 1.  U[1][1]    U[1][1]
    # 2.  U[2][2]    U[1][2]
    # 3.  U[3][3]    U[1][3]
    # 4.  U[1][2]    U[2][2]
    # 5.  U[1][3]    U[2][3]
    # 6.  U[2][3]    U[3][3]
    #

    if ($first_field eq 'ANISOU') {
	$xlat_flag = $xlat_save;

	# parse field

	++$aniso_flag;

	$a_atom_serial_number{$aniso_flag} = substr(($_), 7, 5);
	$aniso_point{($a_atom_serial_number{$aniso_flag} . '|' . $model_flag)}

	  = $aniso_flag;
	$atom_U11{$aniso_flag} = substr(($_), 29, 7);
	$atom_U22{$aniso_flag} = substr(($_), 36, 7);
	$atom_U33{$aniso_flag} = substr(($_), 43, 7);
	$atom_U12{$aniso_flag} = substr(($_), 50, 7);
	$atom_U13{$aniso_flag} = substr(($_), 57, 7);
	$atom_U23{$aniso_flag} = substr(($_), 64, 7);
    }

    #====================================================================
    #  Keyword AUTHOR
    #
    #  Loop over authors as "_audit_author..."

    if ($first_field eq 'AUTHOR') {
	$xlat_flag = $xlat_save;

	# parse record creating an array of authors

	$text = substr(($_), 11, 60);
	$cont = substr(($_), 9, 2);

	{
	    $num_auth = (@authors = split($comma, $text, 9999));
	    if ($num_auth >= 1 && $authors[$num_auth] eq '' &&

	      $comma eq ' ') {
		--$num_auth;
	    }
	}
	for ($i = 1; $i <= $num_auth; ++$i) {
	    {
		$num_a_split = (@a_split = split(' ', $authors[$i], 9999));
		if ($num_a_split >= 1 && $a_split[$num_a_split] eq '' &&

		  ' ' eq ' ') {
		    --$num_a_split;
		}
	    }
	    $authors[$i] = '';
	    if ($num_a_split > 0) {
		$authors[$i] = $a_split[1];
		for ($j = 2; $j <= $num_a_split; ++$j) {
		    $authors[$i] = ($authors[$i] . ' ' . $a_split[$j]);
		}
	    }
	    if ($auth_convtext eq 'yes' ||

	      ($auth_convtext eq 'conditional' && $convtext eq 'yes')) {
		{
		    #
		    # produce a CIF-style name from a PDB name
		    #
		    # begin by applying typesetting codes if any
		    # but always treat "-" and "'"  as breaks for capitalization
		    # in names
		    #
		    {
			$lx_tl = length($authors[$i]);
			$tx_tl = $authors[$i];
			$lostr = '';
			for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
			    $cx_tl = substr($tx_tl, $ix_tl, 1);
			    $cx_tl = substr(($lcaz . $cx_tl),

			      index(($UCAZ . $cx_tl), $cx_tl), 1);
			    $lostr = ($lostr . $cx_tl);
			}
		    }

		    $lstr = length($lostr);
		    $mystr = '';
		    $pchar = ' ';
		    for ($qnsi = 1; $qnsi <= $lstr; ++$qnsi) {
			$mychar = substr($lostr, $qnsi, 1);
			if ($pchar eq ' ' || $pchar eq ',' || $pchar eq '.' ||

			  $pchar eq '-' || $pchar eq "'" || $pchar eq '(' ||

			  $pchar eq '*' || $pchar eq '/') {
			    $mychar = substr(($UCAZ . $mychar),

			      index(($lcaz . $mychar), $mychar), 1);
			}
			if (($mychar ne '*' && $mychar ne "\$" &&

			  $mychar ne '/') || ($mychar eq $pchar)) {	#???
			    $mystr = ($mystr . $mychar);
			}
			if ($pchar eq '/') {
			    if ($mychar eq "\$" || $mychar eq '-') {
				$pchar = $mychar;
			    }
			    # end if( mychar == "$" || mychar == "-" ) 

			    ;
			}
			else {
			    $pchar = $mychar;
			}
			# end if( pchar == "/" )

			;
		    }
		    # end for( qnsi=1; qnsi <= lstr; ++qnsi)
		    #
		    # See if a specific replacement was given
		    #
		    {
			$lx_tu = length($mystr);
			$tx_tu = $mystr;
			$name_temp = '';
			for ($ix_tu = 1; $ix_tu <= $lx_tu; ++$ix_tu) {
			    $cx_tu = substr($tx_tu, $ix_tu, 1);
			    $cx_tu = substr(($UCAZ . $cx_tu),

			      index(($lcaz . $cx_tu), $cx_tu), 1);
			    $name_temp = ($name_temp . $cx_tu);
			}
		    }

		    if ($rep_name{$name_temp} ne '') {
			$mystr = $rep_name{$name_temp};
			#
			# See if there is a comma in place if so we are done
			#
			;
		    }
		    if (index($mystr, $comma) != 0) {
			$ret_val = $mystr;
		    }
		    else {
			$nam_suf = '';
			{
			    $num_namp = (@x_namep = split(' ', $mystr, 9999));
			    if ($num_namp >= 1 && $x_namep[$num_namp] eq '' &&

			      ' ' eq ' ') {
				--$num_namp;
			    }
			}
			if ($num_namp > 1) {
			    {
				$lx_tu = length($x_namep[$num_namp]);
				$tx_tu = $x_namep[$num_namp];
				$xtemp = '';
				for ($ix_tu = 1; $ix_tu <= $lx_tu; ++$ix_tu) {
				    $cx_tu = substr($tx_tu, $ix_tu, 1);
				    $cx_tu = substr(($UCAZ . $cx_tu),

				      index(($lcaz . $cx_tu), $cx_tu), 1);
				    $xtemp = ($xtemp . $cx_tu);
				}
			    }

			    if ($rep_suffix{$xtemp} ne '') {
				if ($junior_on_last eq 'yes') {
				    $x_namep[$num_namp] = $rep_suffix{$xtemp};
				}
				else {
				    $nam_suf = (' ' . $rep_suffix{$xtemp});
				    --$num_namp;
				}
			    }
			    $mystr = $x_namep[1];
			    for ($knamp = 2; $knamp <= $num_namp; ++$knamp) {
				$mystr = ($mystr . ' ' . $x_namep[$knamp]);
			    }
			}
			# end if (num_namp > 1)
			$llname = length($mystr);
			$cc = '';
			for ($kc = $llname - 1; $kc > 1; --$kc) {
			    $cp = $cc;
			    $cc = substr($mystr, $kc, 1);
			    if ($cc eq '.') {	#???
				if ($cp ne ' ') {
				    $mystr = (substr($mystr, $kc + 1,

				      $llname - $kc) . $comma . ' ' .

				      substr($mystr, 1, $kc));
				}
				else {
				    $mystr = (substr($mystr, $kc + 2,

				      $llname - $kc - 1) . $comma . ' ' .

				      substr($mystr, 1, $kc));
				}
				# if (cp != " ")
				$kc = 0;
			    }
			    # end if (cc == ".")

			    ;
			}
			# for (kc=llname-1; kc>1; --kc)
			$mystr = ($mystr . $nam_suf);
			$ret_val = $mystr;
		    }
		    # end if (index(mystr,comma) != 0 )

		    $authors[$i] = $ret_val;
		}
	    }
	}
	$is_blank = $authors[$num_auth];
	if ($is_blank eq '') {
	    --$num_auth;
	}
	if ($num_auth >= 1 && $author_flag eq '1') {	#???
	    printf (("\n\n\n"));
	    printf (("####################\n"));
	    printf (("#                  #\n"));
	    printf (("# AUDIT_AUTHOR     #\n"));
	    printf (("#                  #\n"));
	    printf (("####################\n\n\n"));
	    printf (("loop_ \n"));
	    printf (("_audit_author.name \n"));
	}
	for ($i = 1; $i <= $num_auth; ++$i) {
	    printf "'%s'   \n", $authors[$i];
	}
	if ($num_auth > 0) {
	    ++$author_flag;
	}
    }

    #===========================================================================
    #  Keyword CAVEAT
    #
    # In the 1995 format, a new record, CAVEAT, was added to warn of severe
    # errors in an entry.
    #
    #     caveat_cont     [9-10]
    #     caveat_id       [12-15]
    #     caveat_text     [20-70]       = _struct.title

    if ($first_field eq 'CAVEAT') {
	$xlat_flag = $xlat_save;

	$caveat_cont = substr(($_), 9, 2);
	$caveat_id = substr(($_), 12, 4);
	$caveat_text = substr(($_), 20, 51);

	if ($caveat_cont eq '  ') {
	    $compnd{$compnd_flag++} = 'Warning of Severe Error::';
	}
	$bp = '     ';
	if ($caveat_cont ne '  ') {
	    $bp = '    ';
	}
	$compnd{$compnd_flag++} = ($bp . $caveat_text);
    }

    #===========================================================================
    #   Keyword CISPEP
    #
    #   Introduced with the February 1996 PDB format
    #
    #       cp_sernum              [ 8-10]     
    #       cp_res_name_beg        [12-14]     
    #       cp_chain_id_beg           [16]     
    #       cp_res_seq_num_beg     [18-21]     
    #       cp_icode_beg              [22]     
    #       cp_res_name_end        [26-28]     = _struct_mon_prot.label_comp_id
    #                                            _struct_mon_prot_cis.label_comp_id
    #       cp_chain_id_end           [30]     = _struct_mon_prot.label_asym_id
    #                                            _struct_mon_prot_cis.label_asym_id
    #       cp_res_seq_num_end     [32-35]     = _struct_mon_prot.auth_seq_id
    #                                            _struct_mon_prot_cis.auth_seq_id
    #       cp_icode_end              [36]     append to
    #                                            _struct_mon_prot.auth_seq_id
    #                                            _struct_mon_prot_cis.auth_seq_id
    #       cp_modnum              [44-46]     = _struct_mon_prot.label_model_id
    #                                            _struct_mon_prot_cis.label_model_id
    #       cp_omega               [54-59]     = _struct_mon_prot.omega

    if ($first_field eq 'CISPEP') {
	$xlat_flag = $xlat_save;

	$cp_res_name_end{++$cispep_flag} = substr(($_), 26, 3);
	$cp_chain_id_end{$cispep_flag} = substr(($_), 30, 1);
	$cp_res_seq_num_end{$cispep_flag} = substr(($_), 32, 5);
	$cp_modnum{$cispep_flag} = substr(($_), 44, 3);
	$cp_omega{$cispep_flag} = substr(($_), 54, 6);

	if ($cp_res_name_end{$cispep_flag} eq '   ') {
	    $cp_res_name_end{$cispep_flag} = ' . ';
	}
	if ($cp_chain_id_end{$cispep_flag} eq ' ') {
	    $cp_chain_id_end{$cispep_flag} = $bcid;
	}
	if ($cp_modnum{$cispep_flag} eq '   ') {
	    $cp_modnum{$cispep_flag} = ' . ';
	}
	if ($cp_omega{$cispep_flag} eq '      ') {
	    $cp_omega{$cispep_flag} = '   .   ';
	}
    }

    #==========================================================================
    #  keyword CRYST1
    #
    #

    if ($first_field eq 'CRYST1') {
	$xlat_flag = $xlat_save;
	#
	#  Contains a b c alpha beta gamma SG Z
	#

	# calculate cell volume

	{
	    $ca = cos(substr(($_), 34, 7) * 0.0174532);
	    $cb = cos(substr(($_), 41, 7) * 0.0174532);
	    $cc = cos(substr(($_), 48, 7) * 0.0174532);
	    $cz = (1.0 - ($ca * $ca - $cb * $cb - $cc * $cc) + (2.0 * $ca *

	      $cb * $cc));
	    $vol = (substr(($_), 7, 9) * substr(($_), 16, 9) * substr(($_),

	      25, 9) * (sqrt($cz)));
	    if ($vol - 1 < .01) {
		$warning_list{++$warning_flag} =

		  "#=# CELL: The volume is 1, may be model or NMR, read REMARKs\n";
		++$vol_flag;
	    }
	}
	# localize space group and Z

	{
	    $sg = substr(($_), 56, 11);
	    $Z = substr(($_), 67, 4);
	}

	printf (("\n"));
	printf "_cell.entry_id           %s\n", $head_PDB_code;
	printf "_cell.length_a          %9.3f\n", substr(($_), 7, 9);
	printf "_cell.length_b          %9.3f\n", substr(($_), 16, 9);
	printf "_cell.length_c          %9.3f\n", substr(($_), 25, 9);
	printf "_cell.angle_alpha        %7.2f\n", substr(($_), 34, 7);
	printf "_cell.angle_beta         %7.2f\n", substr(($_), 41, 7);
	printf "_cell.angle_gamma        %7.2f\n", substr(($_), 48, 7);
	printf "_cell.volume         %10.1f \n", $vol;
	printf (("_cell.details                 ? \n"));
	printf "_cell.Z_PDB               %3d \n\n", $Z;
	printf "_symmetry.entry_id                %s \n", $head_PDB_code;
	printf "_symmetry.space_group_name_H-M   '%11s' \n\n", $sg;

	if ($verbose eq 'yes') {
	    printf (("_cell_measurement.temp     ? \n"));
	    printf (("_cell_measurement.theta_min       ? \n"));
	    printf (("_cell_measurement.theta_max       ? \n"));
	    printf (("_cell_measurement.wavelength      ? \n"));
	    printf (("_cell_measurement.pressure        ? \n"));
	    printf (("_cell_measurement.radiation       ? \n"));
	    printf (("_cell_measurement.reflns_used     ? \n\n"));

	    printf (("loop_\n"));
	    printf (("_cell_measurement_refln.index_h \n"));
	    printf (("_cell_measurement_refln.index_k \n"));
	    printf (("_cell_measurement_refln.index_l \n"));
	    printf (("_cell_measurement_refln.theta \n"));
	    printf (("   ?   ?   ?   ? \n"));
	}
    }

    #======================================================================
    #   Keyword COMPND
    #
    # This is considered a common name for the macromolecule
    # in the 1992 format, and a more detailed description with
    # keywords in the 1995 format.  In either case the entire
    # COMPND record is added to the information used for 
    # _struct.title along with the information from TITLE,
    # SOURCE and CAVEAT
    #
    #    record name           [ 1 -  6]     =    "COMPND"
    #    continuation flag     [ 9 - 10]     =    blank for first record
    #    compound              [11 - 70]     =    _struct.title
    #
    #

    if ($first_field eq 'COMPND') {
	$xlat_flag = $xlat_save;
	$compnd_contin = substr(($_), 9, 2);
	if ($compnd_contin eq '  ') {
	    $compnd{$compnd_flag++} = 'Compound::';
	}
	$bp = '     ';
	if ($compnd_contin ne '  ') {
	    $bp = '    ';
	}
	$compnd{$compnd_flag} = ($bp . substr(($_), 11, 60));
	# typeset information, if requested
	if ($convtext eq 'yes') {
	    #
	    # apply PDB typsetting codes if any to a line
	    #
	    {
		$lx_tl = length($compnd{$compnd_flag});
		$tx_tl = $compnd{$compnd_flag};
		$lostr = '';
		for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
		    $cx_tl = substr($tx_tl, $ix_tl, 1);
		    $cx_tl = substr(($lcaz . $cx_tl), index(($UCAZ . $cx_tl),

		      $cx_tl), 1);
		    $lostr = ($lostr . $cx_tl);
		}
	    }

	    $lstr = length($lostr);
	    $mystr = '';
	    $pchar = ' ';
	    for ($qtsi = 1; $qtsi <= $lstr; ++$qtsi) {
		$mychar = substr($lostr, $qtsi, 1);
		if ($pchar eq ' ' || $pchar eq ',' || $pchar eq '.' ||

		  $pchar eq '(' || $pchar eq '*' || $pchar eq '/') {
		    $mychar = substr(($UCAZ . $mychar),

		      index(($lcaz . $mychar), $mychar), 1);
		}
		if (($mychar ne '*' && $mychar ne "\$" && $mychar ne '/') ||

		  ($mychar eq $pchar)) {	#???
		    $mystr = ($mystr . $mychar);
		}
		if ($pchar eq '/') {
		    if ($mychar eq "\$" || $mychar eq '-') {
			$pchar = $mychar;
		    }
		}
		else {
		    $pchar = $mychar;
		}
	    }
	    $ret_val = $mystr;

	    $compnd{$compnd_flag} = $ret_val;
	}

	++$compnd_flag;
    }

    #======================================================================
    #   Keyword CONECT
    #
    #   Origin serial number	= _struct_conn.ptnr1_label_comp_id
    #       			= _struct_conn.ptnr1_label_asym_id
    #       			= _struct_conn.ptnr1_auth_seq_id
    #       			= _struct_conn.ptnr1_label_atom_id
    #       			= _struct_conn.ptnr1_label_alt_id
    #   Target serial numbers	= _struct_conn.ptnr2_label_comp_id
    #       			= _struct_conn.ptnr2_label_asym_id
    #       			= _struct_conn.ptnr2_auth_seq_id
    #       			= _struct_conn.ptnr2_label_atom_id
    #       			= _struct_conn.ptnr2_label_alt_id
    #   Hydrogen bond donor         = _struct_conn.conn_type_id
    #   Hydrogen bond acceptor	= _struct_conn.conn_type_id
    #   Salt bridge excess -ve	= _struct_conn.conn_type_id
    #   Salt bridge excess +ve	= _struct_conn.conn_type_id
    #
    #  _struct_conn.id              = incremental number assigned to each
    #                                 CONECT record
    #  _struct_conn.conn_type_id    = matches generic _struct_conn_type.criteria
    #
    #  all atoms at 1_555 ie no support for -ve targets
    #  No special details included
    #

    if ($first_field eq 'CONECT') {
	$xlat_flag = $xlat_save;
	$connect_save{++$connect_flag} = substr(($_), 1, 61);

	++$conect_flag_2;
    }

    #===========================================================================
    #  Keyword DBREF
    #
    # In the 1995 format, a new record, DBREF, was added to provide
    # "cross-reference links between PDB and the corresponding sequence
    # database entries."  The citations may be to subchains specified
    # by PDB sequence number and insertion code ranges.  
    #
    #     DBREF            [1- 5]
    #     dbref_idcode     [8-11]       = idcode of this entry 
    #     dbref_chainID    [ 13 ]       = _struct_asym.id
    #                                   = _struct_ref.biol_id
    #     dbref_seqBegin  [15-19]       = combines seqBegin and insertBegin
    #                   used to obtain start point in _entity_poly_seq.num
    #                   then mapped to _struct_ref_seq.seq_align_beg
    #     dbref_seqEnd    [21-25]       = combines seqEnd and insertEnd 
    #                   used to obtain start point in _entity_poly_seq.num
    #                   then mapped to _struct_ref_seq.seq_align_end
    #     dbref_database  [27-32]       = _struct_ref.db_name
    #     dbref_dbAccession
    #                     [34-41]       = _struct_ref.db_code
    #     dbref_dbIdCode  [43-54]       = add to _struct_ref.db_code
    #     dbref_dbseqBeg  [56-61]       = _struct_ref_seq.db_align_beg
    #     dbref_dbseqEnd  [63-68]       = _struct_ref_seq.db_align_end
    #
    #     Note:  as of this writing, _struct_ref_seq_dif.db_seq_num is
    #     not in the mmCIF dictionary.
    #
    #     if the database is PDB, columns 61 and 68 contain an insertion code
    #     for other databases, these columns are blank

    if ($first_field eq 'DBREF') {
	$xlat_flag = $xlat_save;
	$dbref_chainID{++$dbref_flag} = substr(($_), 13, 1);
	if ($dbref_chainID{$dbref_flag} eq ' ') {
	    $dbref_chainID{$dbref_flag} = $bcid;
	}
	$dbref_seqBegin{$dbref_flag} = substr(($_), 15, 5);
	$dbref_seqEnd{$dbref_flag} = substr(($_), 21, 5);
	$dbref_database{$dbref_flag} = substr(($_), 27, 6);
	$dbref_dbAccession{$dbref_flag} = substr(($_), 34, 8);
	$dbref_dbIdCode{$dbref_flag} = substr(($_), 43, 12);
	$dbref_dbseqBeg{$dbref_flag} = substr(($_), 56, 6);
	$dbref_dbseqEnd{$dbref_flag} = substr(($_), 63, 6);
	{
	    $numx = (@dblist = split(' ', $dbref_database{$dbref_flag},

	      9999));
	    if ($numx >= 1 && $dblist[$numx] eq '' && ' ' eq ' ') {
		--$numx;
	    }
	}
	$dbref_database{$dbref_flag} = '.';
	if ($numx > 0) {
	    $dbref_database{$dbref_flag} = $dblist[1];
	    for ($j = 2; $j <= $numx; ++$j) {
		$dbref_database{$dbref_flag} = ($dbref_database{$dbref_flag} .

		  '_' . $dblist[$j]);
	    }
	}
	{
	    $numx = (@dblist = split(' ', $dbref_dbAccession{$dbref_flag},

	      9999));
	    if ($numx >= 1 && $dblist[$numx] eq '' && ' ' eq ' ') {
		--$numx;
	    }
	}
	$dbref_dbAccession{$dbref_flag} = '.';
	if ($numx > 0) {
	    $dbref_dbAccession{$dbref_flag} = $dblist[1];
	    for ($j = 2; $j <= $numx; ++$j) {
		$dbref_dbAccession{$dbref_flag} =

		  ($dbref_dbAccession{$dbref_flag} . '_' . $dblist[$j]);
	    }
	}
	{
	    $numx = (@dblist = split(' ', $dbref_dbIdCode{$dbref_flag},

	      9999));
	    if ($numx >= 1 && $dblist[$numx] eq '' && ' ' eq ' ') {
		--$numx;
	    }
	}
	$dbref_dbIdCode{$dbref_flag} = '';
	if ($numx > 0) {
	    $dbref_dbIdCode{$dbref_flag} = (' ' . $dblist[1]);
	    for ($j = 2; $j <= $numx; ++$j) {
		$dbref_dbIdCode{$dbref_flag} = ($dbref_dbIdCode{$dbref_flag} .

		  '_' . $dblist[$j]);
	    }
	}
	if ((' ' . $dbref_dbAccession{$dbref_flag}) eq	#???

	  $dbref_dbIdCode{$dbref_flag}) {
	    $dbref_dbIdCode{$dbref_flag} = '';
	}
    }

    #==========================================================================
    #  keyword END
    #
    #       terminates processing of records, but remainder of file is read
    #
    #

    if ($first_field eq 'END') {
	$xlat_flag = $xlat_save;
	++$end_flag;
    }

    #=============================================================================
    #   Keyword ENDMDL

    if ($first_field eq 'ENDMDL') {
	$xlat_flag = $xlat_save;
	$model_flag = '.';
    }

    #====================================================================
    #   Keyword EXPDTA
    #
    #   expdta [11-70] = _exptl.method
    #

    if ($first_field eq 'EXPDTA') {
	$xlat_flag = $xlat_save;

	# parse field

	$expdta = substr(($_), 11, 60);
	$nmr_flag = index($expdta, 'NMR');
	{
	    $lx_tl = length($expdta);
	    $tx_tl = $expdta;
	    $loexpdta = '';
	    for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
		$cx_tl = substr($tx_tl, $ix_tl, 1);
		$cx_tl = substr(($lcaz . $cx_tl), index(($UCAZ . $cx_tl),

		  $cx_tl), 1);
		$loexpdta = ($loexpdta . $cx_tl);
	    }
	}

	{
	    $num_expdta = (@exp_split = split(' ', $loexpdta, 9999));
	    if ($num_expdta >= 1 && $exp_split[$num_expdta] eq '' &&

	      ' ' eq ' ') {
		--$num_expdta;
	    }
	}
	$loexpdta = '';
	if ($num_expdta > 0) {
	    $loexpdta = $exp_split[1];
	    for ($j = 2; $j <= $num_expdta; ++$j) {
		$loexpdta = ($loexpdta . ' ' . $exp_split[$j]);
	    }
	}
	$expwarn = 'true';
	if ($loexpdta eq 'x-ray diffraction') {
	    $loexpdta = 'single-crystal x-ray diffraction';
	    $expwarn = 'false';
	}
	if ($loexpdta eq 'theoretical model') {
	    $expwarn = 'false';
	}
	if ($expwarn eq 'true') {
	    $warning_list{++$warning_flag} =

	      ('#=# EXPTL: Non-enumerated method: ' . $loexpdta . "\n");
	}

	printf "_exptl.entry_id %s\n", $head_PDB_code;
	printf "_exptl.method  '%-s'\n", $loexpdta;
    }

    #======================================================================
    #   Keyword FORMUL - chemical formula of non-standard groups
    #
    #              component number == _entity.id & _chem_comp.entity_id
    #              het identifier   == _entity_name_common & _chem_comp.id
    #              het_formula_mw   == ignored
    #              het_formula_text == _chem_comp.formula
    #              ??               == _entity_special_details
    #
    #   Information written in ATOM/HETATM keyword

    if ($first_field eq 'FORMUL') {
	$xlat_flag = $xlat_save;

	# parse field

	$formul_het_number{$formul_flag} = substr(($_), 9, 2) + 0;
	$formul_het_site_symbol{$formul_flag} = substr(($_), 13, 3);
	$hetatm_entity = $formul_het_site_symbol{$formul_flag};
	$formul_het_cont_flag{$formul_flag} = substr(($_), 17, 2);
	$hetatm_entity = substr(($_), 13, 3);
	$entity_seq_num{$hetatm_entity} = $formul_het_number{$formul_flag} +

	  0;
	$formul_het_text{$formul_flag} = substr(($_), 20, 51);
	if (substr(($_), 17, 2) eq '  ') {
	    $het_formula{$hetatm_entity} = ("\n; " .

	      $formul_het_text{$formul_flag});
	}
	else {
	    $het_formula{$hetatm_entity} = ($het_formula{$hetatm_entity} .

	      "\n  " . $formul_het_text{$formul_flag});
	}
	$pxxc = '';
	for ($ii = 0; $ii <= 50; ++$ii) {
	    $xxc = substr(($_), 20 + $ii, 1);
	    if (index($UCAZ, $xxc) == 0) {
		$xxc = ' ';
	    }
	    if ($pxxc ne ' ' || $xxc ne ' ') {
		$xhet_formula{$hetatm_entity} = ($xhet_formula{$hetatm_entity}

		  . $xxc);
	    }
	    $pxxc = $xxc;
	}
	$xhet_formula{$hetatm_entity} = ($xhet_formula{$hetatm_entity} . ' ');
	$ent_non_poly_id{$hetatm_entity}++;
	if ($ent_non_poly_id{$hetatm_entity} == 1) {
	    $next_non_poly_id = $ent_non_poly_point{' '};
	    $prev_non_poly_id = ' ';
	    while ($next_non_poly_id ne '') {
		$prev_non_poly_id = $next_non_poly_id;
		$next_non_poly_id = $ent_non_poly_point{$prev_non_poly_id};
	    }
	    $ent_non_poly_point{$prev_non_poly_id} = $hetatm_entity;
	    $ent_non_poly_point{$hetatm_entity} = '';
	    ++$num_non_poly_ents;
	    $ent_non_poly_num{$hetatm_entity} =

	      $formul_het_number{$formul_flag};
	}
	++$formul_flag;

	#  Set up to read addiional entities from ATOM records
	($entity_flag = $formul_flag - 1);
    }

    #=========================================================================
    #   keyword FTNOTE -- footnote to atoms or residues
    #
    #   footnote number  == _atom_sites_footnote.id
    #   footnote text    == _atom_sites_footnote.text

    if ($first_field eq 'FTNOTE') {
	$xlat_flag = $xlat_save;

	$X = substr(($_), 10, 1);

	if ($ftnote_flag_old + 0 == -1) {
	    $ft_save{++$foot_flag} = "\nloop_\n";
	    $ft_save{++$foot_flag} = "_atom_sites_footnote.id  \n";
	    $ft_save{++$foot_flag} = "_atom_sites_footnote.text \n";
	}
	$ftnote_num = substr(($_), 8, 3) + 0;
	$ftnote_text = substr(($_), 12, 59);

	if ($ftnote_num != $ftnote_flag_old) {	#???
	    if ($ftnote_flag_old + 0 != -1) {
		$ft_save{++$foot_flag} = "; \n";
	    }
	    $ft_save{++$foot_flag} = ('  ' . $ftnote_num . "\n");
	    $ft_save{++$foot_flag} = ('; ' . $ftnote_text . "\n");
	    ++$ftnote_flag{$ftnote_num};
	}
	else {
	    $ft_save{++$foot_flag} = ('  ' . $ftnote_text . "\n");
	    ++$ftnote_flag{$ftnote_num};
	}
	$ftnote_flag_old = $ftnote_num;
    }

    #====================================================================
    #  Keyword HEADER
    #
    #  This is a good place to place the _struct_biol data items. Templates
    #  are given but no information has been parsed excluding
    #  _special_details.
    #
    #       head_funct_class  [11-50]   ==  _struct_biol.details
    #       head_dep_date     [51-59]   ==  _database_PDB_rev.date_original
    #                                       _audit.creation_date
    #       head_PDB_code     [63-66]   ==  _database_2.database_code
    #                                       _struct_biol.id
    #                                       _audit.revision_id

    if ($first_field eq 'HEADER') {
	$xlat_flag = $xlat_save;

	$head_funct_class = substr(($_), 11, 40);
	$head_dep_date = substr(($_), 51, 9);
	$head_PDB_code = substr(($_), 63, 4);
	if ($head_PDB_code eq '    ') {
	    $head_PDB_code = '.';
	}
	#
	#       Output the PDB code immediately as the data block name
	#
	printf "data_%4s\n\n", $head_PDB_code;
	printf "_entry.id        %4s\n\n", $head_PDB_code;
	#
	#       save the header id as a possible audit.revision_id
	#
	$aud_rev_id = $head_PDB_code;
    }

    #=======================================================================
    #   Keyword HELIX
    #
    #  8 - 10 helix_no.             == (not used)
    # 12 - 14 helix_id              == _struct_conf.id
    # 16 - 18 helix_res_name_beg    == _struct_conf.beg_label_comp_id
    # 20      helix_chain_id_beg    == _struct_conf.beg_label_asym_id
    # 22 - 26 helix_res_seq_beg     == _struct_conf.beg_auth_seq_id
    # 28 - 30 helix_res_name_end    == _struct_conf.end_label_comp_id
    # 32      helix_chain_id_end    == _struct_conf.end_label_asym_id
    # 34 - 38 helix_res_seq_end     == _struct_conf.end_auth_seq_id
    # 39 - 40 helix_class           == _struct_conf.conf_type_id
    # 41 - 70 helix_comment         == _struct_conf.details
    #
    # note helix classes 9 and 10 as defined by the PDB do not have CIF 
    # definitions
    #
    #

    if ($first_field eq 'HELIX') {
	$xlat_flag = $xlat_save;

	$helix_no{$ss_flag} = substr(($_), 8, 3);
	$helix_id{$ss_flag} = substr(($_), 12, 3);
	$helix_res_name_beg{$ss_flag} = substr(($_), 16, 3);
	$helix_chain_id_beg{$ss_flag} = substr(($_), 20, 1);
	$helix_res_seq_beg{$ss_flag} = substr(($_), 22, 5);
	$helix_res_name_end{$ss_flag} = substr(($_), 28, 3);
	$helix_chain_id_end{$ss_flag} = substr(($_), 32, 1);
	$helix_res_seq_end{$ss_flag} = substr(($_), 34, 5);
	$helix_class{$ss_flag} = substr(($_), 39, 2);
	$helix_comment{$ss_flag} = substr(($_), 41, 30);
	if ($helix_comment{$ss_flag} eq '                              ' ||

	  $helix_comment{$ss_flag} eq '') {
	    $helix_comment{$ss_flag} = '  .  ';
	    if ($helix_class{$ss_flag} + 0 == 1) {
		$helix_comment{$ss_flag} = 'RIGHT-HANDED ALPHA HELIX';
	    }
	    if ($helix_class{$ss_flag} + 0 == 2) {
		$helix_comment{$ss_flag} = 'RIGHT-HANDED OMEGA HELIX';
	    }
	    if ($helix_class{$ss_flag} + 0 == 3) {
		$helix_comment{$ss_flag} = 'RIGHT-HANDED PI HELIX';
	    }
	    if ($helix_class{$ss_flag} + 0 == 4) {
		$helix_comment{$ss_flag} = 'RIGHT-HANDED GAMMA HELIX';
	    }
	    if ($helix_class{$ss_flag} + 0 == 5) {
		$helix_comment{$ss_flag} = 'RIGHT-HANDED 3/10 HELIX';
	    }
	    if ($helix_class{$ss_flag} + 0 == 6) {
		$helix_comment{$ss_flag} = 'LEFT-HANDED ALPHA HELIX';
	    }
	    if ($helix_class{$ss_flag} + 0 == 7) {
		$helix_comment{$ss_flag} = 'LEFT-HANDED OMEGA HELIX';
	    }
	    if ($helix_class{$ss_flag} + 0 == 8) {
		$helix_comment{$ss_flag} = 'LEFT-HANDED GAMMA HELIX';
	    }
	    if ($helix_class{$ss_flag} + 0 == 9) {
		$helix_comment{$ss_flag} = '2/7 RIBBON/HELIX';
	    }
	    if ($helix_class{$ss_flag} + 0 == 10) {
		$helix_comment{$ss_flag} = 'POLYPROLINE';
	    }
	}

	# strip blanks from id
	{
	    $num_x = (@xxx = split(' ', $helix_id{$ss_flag}, 9999));
	    if ($num_x >= 1 && $xxx[$num_x] eq '' && ' ' eq ' ') {
		--$num_x;
	    }
	}
	$helix_id{$ss_flag} = '';
	if ($num_x == 1) {
	    $helix_id{$ss_flag} = $xxx[1];
	}
	if ($num_x == 2) {
	    $helix_id{$ss_flag} = ($xxx[1] . '_' . $xxx[2]);
	    # provide default conditions

	    ;
	}
	if ($helix_chain_id_beg{$ss_flag} eq ' ') {
	    $helix_chain_id_beg{$ss_flag} = $bcid;
	}
	if ($helix_chain_id_end{$ss_flag} eq ' ') {
	    $helix_chain_id_end{$ss_flag} = $bcid;

	    # give real names to helix classes

	    ;
	}
	$h_class_suffix = '_P';
	{
	    $num_x = (@xxx = split(' ',

	      ($helix_res_name_end{$ss_flag} . ' ' .

	      $helix_res_name_beg{$ss_flag}), 9999));
	    if ($num_x >= 1 && $xxx[$num_x] eq '' && ' ' eq ' ') {
		--$num_x;
	    }
	}
	foreach $i ($[ .. $#na_list) {
	    if ($na_list[$i] eq $xxx[1] || $na_list[$i] eq $xxx[2]) {	#???	#???
		$h_class_suffix = '_N';
	    }
	}

	if ($helix_class{$ss_flag} eq ' 1') {
	    $helix_class{$ss_flag} = ('HELX_RH_AL' . $h_class_suffix);
	}
	if ($helix_class{$ss_flag} eq ' 2') {
	    $helix_class{$ss_flag} = ('HELX_RH_OM' . $h_class_suffix);
	}
	if ($helix_class{$ss_flag} eq ' 3') {
	    $helix_class{$ss_flag} = ('HELX_RH_PI' . $h_class_suffix);
	}
	if ($helix_class{$ss_flag} eq ' 4') {
	    $helix_class{$ss_flag} = ('HELX_RH_GA' . $h_class_suffix);
	}
	if ($helix_class{$ss_flag} eq ' 5') {
	    $helix_class{$ss_flag} = ('HELX_RH_3T' . $h_class_suffix);
	}
	if ($helix_class{$ss_flag} eq ' 6') {
	    $helix_class{$ss_flag} = ('HELX_LH_AL' . $h_class_suffix);
	}
	if ($helix_class{$ss_flag} eq ' 7') {
	    $helix_class{$ss_flag} = ('HELX_LH_OM' . $h_class_suffix);
	}
	if ($helix_class{$ss_flag} eq ' 8') {
	    $helix_class{$ss_flag} = ('HELX_LH_GA' . $h_class_suffix);
	}
	if ($helix_class{$ss_flag} eq ' 9') {
	    $helix_class{$ss_flag} = ('HELX_27' . $h_class_suffix);
	}
	if ($helix_class{$ss_flag} eq '10') {
	    $helix_class{$ss_flag} = ('HELX_PP' . $h_class_suffix);
	}
	if (substr($helix_class{$ss_flag}, 1, 1) ne 'H') {
	    $warning_list{++$warning_flag} =

	      ('#=# STRUCT_CONF: Unrecognized helix classification: ' .

	      $helix_class{$ss_flag} . "\n");
	    $helix_class{$ss_flag} = '.';
	}

	++$ss_flag;
	++$helix_flag;
    }

    #===================================================================
    #   Keyword HET
    #
    #     het_site_symbol      [ 8 - 10] == to link to _entity.id from FORMUL
    #     het_site_chain       [13]      == ????
    #     het_site_seqNum      [14 - 17] == sequence no, or -999 if more than 15
    #     het_site_iCode       [18]      == append to seq no
    #     het_atoms_number     [21 - 25] == THIS IS THE NUMBER OF HETATM LINES
    #                                       NOT A COUNT OF ATOMS
    #                                       By careful processing of the HETATM
    #                                       information, paying attention to
    #                                       occupancies, this number could
    #                                       be related to
    #                                         _chem_comp.number_atoms_all or
    #                                         _chem_comp.number_atoms_nh
    #                                       but we do not attempt this
    #     het_site_text        [31-70]  == _chem_comp.details
    #
    #

    if ($first_field eq 'HET') {
	$xlat_flag = $xlat_save;

	# parse field

	$het_site_symbol{$het_flag} = substr(($_), 8, 3);
	$het_site_chain{$het_flag} = substr(($_), 13, 1);
	$het_site_residue{$het_flag} = substr(($_), 14, 5);
	$het_atoms_number{$het_flag} = substr(($_), 21, 4);
	$het_site_text{$het_site_symbol{$het_flag}} = substr(($_), 31, 40);

	if ($het_site_chain{$het_flag} eq ' ') {
	    $het_site_chain{$het_flag} = '.';
	}

	++$het_flag;
    }

    #===================================================================
    #   Keyword HETNAM
    #
    #     hetnam_cont          [ 9 - 10] == continuation flag
    #     hetnam_symbol        [12 - 14] == to link to entity_id from FORMUL
    #     hetnam_text          [16 - 70] == text of chemical name
    #                                         _chem_comp.name and
    #                                         _entity_name_com.name
    #

    if ($first_field eq 'HETNAM') {
	$xlat_flag = $xlat_save;

	# parse field

	$hetnam_cont{$hetnam_flag} = substr(($_), 9, 2);
	$hetnam_symbol{$hetnam_flag} = substr(($_), 12, 3);
	$hetnam_text{$hetnam_flag} = substr(($_), 16, 55);
	if ($hetnam_cont{$hetnam_flag} eq '  ') {
	    $het_site_name{$hetnam_symbol{$hetnam_flag}} = (' ' . substr(($_),

	      16, 55));
	}
	else {
	    $het_site_name{$hetnam_symbol{$hetnam_flag}} =

	      ($het_site_name{$hetnam_symbol{$hetnam_flag}} . "\n  " .

	      substr(($_), 16, 55));
	}

	++$hetnam_flag;
    }

    #===================================================================
    #   Keyword HETSYN
    #
    #     hetsyn_cont          [ 9 - 10] == continuation flag
    #     hetsyn_symbol        [12 - 14] == to link to entity_id from FORMUL
    #     hetsyn_text          [16 - 70] == text of chemical name
    #                                         _entity_name_com.name
    #

    if ($first_field eq 'HETSYN') {
	$xlat_flag = $xlat_save;

	# parse field

	$hetsyn_cont{$hetsyn_flag} = substr(($_), 9, 2);
	$hetsyn_symbol{$hetsyn_flag} = substr(($_), 12, 3);
	$hetsyn_text{$hetsyn_flag} = substr(($_), 16, 55);
	if ($hetsyn_cont{$hetsyn_flag} eq '  ') {
	    $het_site_syn{$hetsyn_symbol{$hetsyn_flag}} = (substr(($_), 16,

	      55));
	}
	else {
	    $het_site_syn{$hetsyn_symbol{$hetsyn_flag}} =

	      ($het_site_syn{$hetsyn_symbol{$hetsyn_flag}} . substr(($_), 16,

	      55));
	}

	++$hetsyn_flag;
    }

    #===========================================================================
    #   Keyword HYDBND
    #
    #   Introduced with the February 1996 PDB format
    #
    #   There is no way to define the hydrogen atom of a hydrogen bond in mmCIF
    #   we treat it as a partner in a hydrogen bond with role "hydrogen"
    #
    #
    #       hb_atom_beg            [13-16]     = _struct_conn.ptnr1_label_atom_id
    #       hb_alt_loc_beg            [17]     = _struct_conn.ptnr1_label_alt_id
    #       hb_res_name_beg        [18-20]     = _struct_conn.ptnr1_label_comp_id
    #       hb_chain_id_beg           [22]     = _struct_conn.ptnr1_label_asym_id
    #       hb_res_seq_num_beg     [24-27]     = _struct_conn.ptnr1_auth_seq_id
    #       hb_icode_beg              [28]     append to
    #                                            _struct_conn.ptnr1_auth_seq_id
    #       hb_name_ha             [30-33]     = _struct_conn.ptnr1_label_atom_id
    #       hb_alt_loc_ha             [34]     = _struct_conn.ptnr1_label_alt_id
    #       hb_chain_id_ha            [36]     = _struct_conn.ptnr1_label_asym_id
    #       hb_res_seq_num_ha      [37-41]     = _struct_conn.ptnr1_auth_seq_id
    #       hb_icode_ha               [42]     append to
    #                                            _struct_conn.ptnr1_auth_seq_id
    #       hb_atom_end            [44-47]     = _struct_conn.ptnr2_label_atom_id
    #       hb_alt_loc_end            [48]     = _struct_conn.ptnr2_label_alt_id
    #       hb_res_name_end        [49-51]     = _struct_conn.ptnr2_label_comp_id
    #       hb_chain_id_end           [53]     = _struct_conn.ptnr2_label_asym_id
    #       hb_res_seq_num_end     [55-58]     = _struct_conn.ptnr2_auth_seq_id
    #       hb_icode_end              [59]     append to
    #                                            _struct_conn.ptnr2_auth_seq_id
    #       hb_symop1              [60-65]     = _struct_conn.ptnr1_symmetry
    #       hb_symop2              [67-72]     = _struct_conn.ptnr2_symmetry
    #

    if ($first_field eq 'HYDBND') {
	$xlat_flag = $xlat_save;

	$hb_atom_beg{++$hydbnd_flag} = substr(($_), 13, 4);
	$hb_alt_loc_beg{$hydbnd_flag} = substr(($_), 17, 1);
	$hb_res_name_beg{$hydbnd_flag} = substr(($_), 18, 3);
	$hb_chain_id_beg{$hydbnd_flag} = substr(($_), 22, 1);
	$hb_res_seq_num_beg{$hydbnd_flag} = substr(($_), 24, 5);
	$hb_atom_ha{$hydbnd_flag} = substr(($_), 30, 4);
	$hb_alt_loc_ha{$hydbnd_flag} = substr(($_), 34, 1);
	$hb_chain_id_ha{$hydbnd_flag} = substr(($_), 36, 1);
	$hb_res_seq_num_ha{$hydbnd_flag} = substr(($_), 38, 5);
	$hb_atom_end{$hydbnd_flag} = substr(($_), 44, 4);
	$hb_alt_loc_end{$hydbnd_flag} = substr(($_), 48, 1);
	$hb_res_name_end{$hydbnd_flag} = substr(($_), 49, 3);
	$hb_chain_id_end{$hydbnd_flag} = substr(($_), 53, 1);
	$hb_res_seq_num_end{$hydbnd_flag} = substr(($_), 55, 5);
	$hb_symm_1{$hydbnd_flag} = substr(($_), 60, 6);
	$hb_symm_2{$hydbnd_flag} = substr(($_), 67, 6);

	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = $hb_atom_beg{$hydbnd_flag};
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $temp_name = $ret_val;
	}
	$hb_atom_beg{$hydbnd_flag} = $temp_name;
	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = $hb_atom_end{$hydbnd_flag};
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $temp_name = $ret_val;
	}
	$hb_atom_end{$hydbnd_flag} = $temp_name;
	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = $hb_atom_ha{$hydbnd_flag};
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $temp_name = $ret_val;
	}
	$hb_atom_ha{$hydbnd_flag} = $temp_name;
	if ($hb_alt_loc_beg{$hydbnd_flag} eq ' ') {
	    $hb_alt_loc_beg{$hydbnd_flag} = '.';
	}
	if ($hb_alt_loc_end{$hydbnd_flag} eq ' ') {
	    $hb_alt_loc_end{$hydbnd_flag} = '.';
	}
	if ($hb_alt_loc_ha{$hydbnd_flag} eq ' ') {
	    $hb_alt_loc_ha{$hydbnd_flag} = '.';
	}
	if ($hb_res_name_beg{$hydbnd_flag} eq '   ') {
	    $hb_res_name_beg{$hydbnd_flag} = ' . ';
	}
	if ($hb_res_name_end{$hydbnd_flag} eq '   ') {
	    $hb_res_name_end{$hydbnd_flag} = ' . ';
	}
	if ($hb_res_seq_num_beg{$hydbnd_flag} eq '     ') {
	    $hb_res_seq_num_beg{$hydbnd_flag} = '  .  ';
	}
	if ($hb_res_seq_num_end{$hydbnd_flag} eq '     ') {
	    $hb_res_seq_num_end{$hydbnd_flag} = '  .  ';
	}
	if ($hb_res_seq_num_ha{$hydbnd_flag} eq '     ') {
	    $hb_res_seq_num_ha{$hydbnd_flag} = '  .  ';
	}
	if ($hb_chain_id_beg{$hydbnd_flag} eq ' ') {
	    $hb_chain_id_beg{$hydbnd_flag} = $bcid;
	}
	if ($hb_chain_id_end{$hydbnd_flag} eq ' ') {
	    $hb_chain_id_end{$hydbnd_flag} = $bcid;
	}
	if ($hb_chain_id_ha{$hydbnd_flag} eq ' ') {
	    $hb_chain_id_ha{$hydbnd_flag} = $bcid;
	}
	if ($hb_symm_1{$hydbnd_flag} eq '      ') {
	    $hb_symm_1{$hydbnd_flag} = '   .   ';
	}
	else {
	    $hb_symm_1{$hydbnd_flag} = (substr($hb_symm_1{$hydbnd_flag}, 1,

	      3) . '_' . substr($hb_symm_1{$hydbnd_flag}, 4, 3));
	}
	if ($hb_symm_2{$hydbnd_flag} eq '      ') {
	    $hb_symm_2{$hydbnd_flag} = '   .   ';
	}
	else {
	    $hb_symm_2{$hydbnd_flag} = (substr($hb_symm_2{$hydbnd_flag}, 1,

	      3) . '_' . substr($hb_symm_2{$hydbnd_flag}, 4, 3));
	}
    }

    #==================================================================
    #  Keyword JRNL
    #
    # As defined by the PDB, this is the primary citation that matches the
    # given coordinate set. It is written before the REMARK 2 record
    #
    #  "primary"                     =   _citation.id
    #                                =   _citation_author.citation_id
    #  "yes"/"no"                    =   _citation.coordinate_linkage
    #  jrnl_rec_type   [13-16]       =
    #  jrnl_cont       [17-18]       =
    #
    #AUTH
    #  jrnl_auth       [20-70]       =   _citation_author.name
    #
    #TITL
    #  jrnl_titl       [20-70]       =   _citation.title
    #
    #REF
    #  jrnl_ref_jour   [20-47]       =   _citation.journal_abbrev (this is not
    #                                        always abbreviated but it will do)
    #  jrnl_ref_vol    [53-55]       =   _citation.journal_volume
    #      "?"                       =   _citation.journal_issue
    #  jrnl_ref_page   [57-61]       =   _citation.page_first
    #     "?"                        =   _citation.page_last
    #  jrnl_ref_year   [63-66]       =   _citation.year
    #
    #PUBL
    #  jrnl_pub_pub    [20-70]       =  _citation.book_publisher
    #
    #REFN
    #  jrnl_astm       [25-30]       =  _citation.journal_id_ASTM
    #  jrnl_country    [33-34]       =  _citation.country
    #  jrnl_isbn       [41-65]       =  _citation.journal_id_ISSN or
    #                                =  _citation.book_id_ISBN
    #     "?"                        =  _citation.abstract
    #
    #                                =  _citation.details

    if ($first_field eq 'JRNL') {
	$xlat_flag = $xlat_save;
	$flush_refs = 1;

	$jrnl_rec_type = substr(($_), 13, 4);
	$jrnl_cont = substr(($_), 17, 2);
	$jrnl_title = substr(($_), 20, 51);
	$text = substr(($_), 20, 51);
	$jrnl_auth = $text;
	if ($convtext eq 'yes') {
	    #
	    # apply PDB typsetting codes if any to a line
	    #
	    {
		$lx_tl = length($text);
		$tx_tl = $text;
		$lostr = '';
		for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
		    $cx_tl = substr($tx_tl, $ix_tl, 1);
		    $cx_tl = substr(($lcaz . $cx_tl), index(($UCAZ . $cx_tl),

		      $cx_tl), 1);
		    $lostr = ($lostr . $cx_tl);
		}
	    }

	    $lstr = length($lostr);
	    $mystr = '';
	    $pchar = ' ';
	    for ($qtsi = 1; $qtsi <= $lstr; ++$qtsi) {
		$mychar = substr($lostr, $qtsi, 1);
		if ($pchar eq ' ' || $pchar eq ',' || $pchar eq '.' ||

		  $pchar eq '(' || $pchar eq '*' || $pchar eq '/') {
		    $mychar = substr(($UCAZ . $mychar),

		      index(($lcaz . $mychar), $mychar), 1);
		}
		if (($mychar ne '*' && $mychar ne "\$" && $mychar ne '/') ||

		  ($mychar eq $pchar)) {	#???
		    $mystr = ($mystr . $mychar);
		}
		if ($pchar eq '/') {
		    if ($mychar eq "\$" || $mychar eq '-') {
			$pchar = $mychar;
		    }
		}
		else {
		    $pchar = $mychar;
		}
	    }
	    $ret_val = $mystr;

	    $text = $ret_val;
	}

	$jrnl_title = $text;
	$cit_flag = 1;
	$primary = 'yes';
	$cit_refNum{$cit_flag} = 'primary';

	if ($jrnl_rec_type eq 'TITL' && $jrnl_cont eq '  ') {
	    $cit_title_1{$cit_flag} = $jrnl_title;
	    $cit_title_2{$cit_flag} = '';
	}
	if ($jrnl_rec_type eq 'TITL' && $jrnl_cont ne '  ') {
	    if ($cit_title_2{$cit_flag} eq '') {
		$cit_title_2{$cit_flag} = $jrnl_title;
	    }
	    else {
		$cit_title_2{$cit_flag} = ($cit_title_2{$cit_flag} . "\n  " .

		  $jrnl_title);
	    }
	}

	if ($jrnl_rec_type eq 'AUTH' && $jrnl_cont eq '  ') {
	    $cit_auth_1{$cit_flag} = $jrnl_auth;
	    $cit_auth_2{$cit_flag} = '';
	}

	if ($jrnl_rec_type eq 'AUTH' && $jrnl_cont ne '  ') {
	    $cit_auth_2{$cit_flag} = ($cit_auth_2{$cit_flag} . $jrnl_auth);
	}
	if ($jrnl_rec_type eq 'REF ' && $jrnl_cont eq '  ') {
	    $jour_1{$cit_flag} = substr(($_), 20, 28);
	    $jour_2{$cit_flag} = '';
	    $volu{$cit_flag} = substr(($_), 52, 4);
	    $page{$cit_flag} = substr(($_), 57, 5);
	    $year{$cit_flag} = substr(($_), 63, 4);
	    $jrnl_pub_pub_1{$cit_flag} = '?';
	}
	if ($jrnl_rec_type eq 'REF ' && $jrnl_cont ne '  ') {
	    if ($jour_2{$cit_flag} eq '') {
		$jour_2{$cit_flag} = substr(($_), 20, 28);
	    }
	    else {
		$jour_2{$cit_flag} = ($jour_2{$cit_flag} . "\n  " .

		  substr(($_), 20, 28));
	    }
	}
	if ($jrnl_rec_type eq 'PUBL' && $jrnl_cont eq '  ') {
	    $jrnl_pub_pub{$cit_flag} = substr(($_), 20, 51);
	    $jour_1{$cit_flag} = '?';
	    if ($volu{$cit_flag} eq '' || $volu{$cit_flag} eq '    ') {
		$volu{$cit_flag} = '?';
	    }
	    $page{$cit_flag} = '?';
	    $year{$cit_flag} = '?';
	}
	if ($jrnl_rec_type eq 'PUBL' && $jrnl_cont eq '  ') {
	    $jrnl_pub_pub_1{$cit_flag} = substr(($_), 20, 51);
	    $jrnl_pub_pub_2{$cit_flag} = '';
	}
	if ($jrnl_rec_type eq 'PUBL' && $jrnl_cont ne '  ') {
	    if ($jrnl_pub_pub_2{$cit_flag} eq '') {
		$jrnl_pub_pub_2{$cit_flag} = substr(($_), 20, 51);
	    }
	    else {
		$jrnl_pub_pub_2{$cit_flag} = ($jrnl_pub_pub_2{$cit_flag} .

		  "\n  " . substr(($_), 20, 51));
	    }
	}
	if ($jrnl_rec_type eq 'REFN') {
	    $astm{$cit_flag} = substr(($_), 25, 6);
	    $country{$cit_flag} = substr(($_), 33, 2);
	    if ($country{$cit_flag} eq '  ') {
		$country{$cit_flag} = '?';
	    }
	    $issn_isbn{$cit_flag} = substr(($_), 36, 4);
	    if ($issn_isbn{$cit_flag} eq '    ' &&

	      substr($jour_1{$cit_flag}, 1, 9) ne 'TO BE PUB') {
		if ($jrnl_pub_pub_1{$cit_flag} ne '?') {
		    $issn_isbn{$cit_flag} = 'ISBN';
		}
		if ($volu{$cit_flag} eq '' || $volu{$cit_flag} eq '?' ||

		  $volu{$cit_flag} eq '    ') {
		    $issn_isbn{$cit_flag} = 'ISBN';
		}
	    }
	    if ($issn_isbn{$cit_flag} ne 'ISBN') {
		$isbn{$cit_flag} = '?';
		$issn{$cit_flag} = substr(($_), 41, 25);
	    }
	    else {
		$issn{$cit_flag} = '?';
		$isbn{$cit_flag} = substr(($_), 41, 25);
	    }
	    $csd{$cit_flag} = substr(($_), 67, 4);
	    if ($csd{$cit_flag} eq '    ') {
		$csd{$cit_flag} = '?';
	    }
	}
    }

    #====================================================================
    #  Keyword KEYWRD
    #
    #       keywrd_list [11-70]     ==  _struct_keywords.text
    #       from HEADER:
    #         head_PDB_code         ==  _struct_keywords.entry_id

    if ($first_field eq 'KEYWRD' || $first_field eq 'KEYWDS') {
	$xlat_flag = $xlat_save;

	$keywrd_list = substr(($_), 11, 60);
	#
	if ($keywrd_flag == 0) {
	    $keywrd_tail = '';
	    $key_save{++$keywrd_flag} = "\n\n\n";
	    $key_save{++$keywrd_flag} = "#############################\n";
	    $key_save{++$keywrd_flag} = "#                           #\n";
	    $key_save{++$keywrd_flag} = "# STRUCT_KEYWORDS           #\n";
	    $key_save{++$keywrd_flag} = "#                           #\n";
	    $key_save{++$keywrd_flag} = "#############################\n\n";
	    $key_save{++$keywrd_flag} = "loop_\n";
	    $key_save{++$keywrd_flag} = "_struct_keywords.entry_id\n";
	    $key_save{++$keywrd_flag} = "_struct_keywords.text\n";
	}
	if ($key_tail ne '') {
	    --$keywrd_flag;
	    $keywrd_list = ($key_tail . $keywrd_list);
	}
	{
	    $num_keys = (@keyslist = split($comma, $keywrd_list, 9999));
	    if ($num_keys >= 1 && $keyslist[$num_keys] eq '' &&

	      $comma eq ' ') {
		--$num_keys;
	    }
	}
	for ($ii = 1; $ii <= $num_keys; ++$ii) {
	    {
		$num_k_split = (@keys_split = split(' ', $keyslist[$ii],

		  9999));
		if ($num_k_split >= 1 && $keys_split[$num_k_split] eq '' &&

		  ' ' eq ' ') {
		    --$num_k_split;
		}
	    }
	    $keyslist[$ii] = '';
	    if ($num_k_split > 0) {
		$keyslist[$ii] = $keys_split[1];
		for ($j = 2; $j <= $num_k_split; ++$j) {
		    $keyslist[$ii] = ($keyslist[$ii] . ' ' . $keys_split[$j]);
		}
		$key_save{++$keywrd_flag} = sprintf("%4s  '%s'\n",

		  $head_PDB_code, $keyslist[$ii]);
	    }
	    $key_tail = '';
	    if ($num_keys > 0) {
		$key_tail = $keyslist[$num_keys];
	    }
	}
    }

    #===========================================================================
    #   Keyword LINK
    #
    #   Introduced with the February 1996 PDB format
    #
    #       lk_atom_beg            [13-16]     = _struct_conn.ptnr1_label_atom_id
    #       lk_alt_loc_beg            [17]     = _struct_conn.ptnr1_label_alt_id
    #       lk_res_name_beg        [18-20]     = _struct_conn.ptnr1_label_comp_id
    #       lk_chain_id_beg           [22]     = _struct_conn.ptnr1_label_asym_id
    #       lk_res_seq_num_beg     [23-26]     = _struct_conn.ptnr1_auth_seq_id
    #       lk_icode_beg              [27]     append to
    #                                            _struct_conn.ptnr1_auth_seq_id
    #       lk_atom_end            [43-46]     = _struct_conn.ptnr2_label_atom_id
    #       lk_alt_loc_end            [47]     = _struct_conn.ptnr2_label_alt_id
    #       lk_res_name_end        [48-50]     = _struct_conn.ptnr2_label_comp_id
    #       lk_chain_id_end           [52]     = _struct_conn.ptnr2_label_asym_id
    #       lk_res_seq_num_end     [53-56]     = _struct_conn.ptnr2_auth_seq_id
    #       lk_icode_end              [57]     append to
    #                                            _struct_conn.ptnr2_auth_seq_id
    #       lk_symop1              [60-65]     = _struct_conn.ptnr1_symmetry
    #       lk_symop2              [67-72]     = _struct_conn.ptnr2_symmetry
    #

    if ($first_field eq 'LINK') {
	$xlat_flag = $xlat_save;

	$lk_atom_beg{++$link_flag} = substr(($_), 13, 4);
	$lk_alt_loc_beg{$link_flag} = substr(($_), 17, 1);
	$lk_res_name_beg{$link_flag} = substr(($_), 18, 3);
	$lk_chain_id_beg{$link_flag} = substr(($_), 22, 1);
	$lk_res_seq_num_beg{$link_flag} = substr(($_), 23, 5);
	$lk_atom_end{$link_flag} = substr(($_), 43, 4);
	$lk_alt_loc_end{$link_flag} = substr(($_), 47, 1);
	$lk_res_name_end{$link_flag} = substr(($_), 48, 3);
	$lk_chain_id_end{$link_flag} = substr(($_), 52, 1);
	$lk_res_seq_num_end{$link_flag} = substr(($_), 53, 5);
	$lk_symm_1{$link_flag} = substr(($_), 60, 6);
	$lk_symm_2{$link_flag} = substr(($_), 67, 6);

	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = $lk_atom_beg{$link_flag};
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $temp_name = $ret_val;
	}
	$lk_atom_beg{$link_flag} = $temp_name;
	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = $lk_atom_end{$link_flag};
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $temp_name = $ret_val;
	}
	$lk_atom_end{$link_flag} = $temp_name;
	if ($lk_alt_loc_beg{$link_flag} eq ' ') {
	    $lk_alt_loc_beg{$link_flag} = '.';
	}
	if ($lk_alt_loc_end{$link_flag} eq ' ') {
	    $lk_alt_loc_end{$link_flag} = '.';
	}
	if ($lk_res_name_beg{$link_flag} eq '   ') {
	    $lk_res_name_beg{$link_flag} = ' . ';
	}
	if ($lk_res_name_end{$link_flag} eq '   ') {
	    $lk_res_name_end{$link_flag} = ' . ';
	}
	if ($lk_chain_id_beg{$link_flag} eq ' ') {
	    $lk_chain_id_beg{$link_flag} = $bcid;
	}
	if ($lk_chain_id_end{$link_flag} eq ' ') {
	    $lk_chain_id_end{$link_flag} = $bcid;
	}
	if ($lk_symm_1{$link_flag} eq '      ') {
	    $lk_symm_1{$link_flag} = '   .   ';
	}
	else {
	    $lk_symm_1{$link_flag} = (substr($lk_symm_1{$link_flag}, 1,

	      3) . '_' . substr($lk_symm_1{$link_flag}, 4, 3));
	}
	if ($lk_symm_2{$link_flag} eq '      ') {
	    $lk_symm_2{$link_flag} = '   .   ';
	}
	else {
	    $lk_symm_2{$link_flag} = (substr($lk_symm_2{$link_flag}, 1,

	      3) . '_' . substr($lk_symm_2{$link_flag}, 4, 3));
	}
    }

    #=============================================================================
    #   Keyword MASTER
    #
    #   (used in END statement)

    if ($first_field eq 'MASTER') {
	$xlat_flag = $xlat_save;

	# parse totals

	$total_remark = substr(($_), 11, 5);
	$total_ftnote = substr(($_), 16, 5);
	$total_het = substr(($_), 21, 5);
	$total_helix = substr(($_), 26, 5);
	$total_sheet = substr(($_), 31, 5);
	$total_turn = substr(($_), 36, 5);
	$total_site = substr(($_), 41, 5);
	$total_o_s_m = substr(($_), 46, 5);
	$total_a_h = substr(($_), 51, 5);
	$total_ter = substr(($_), 56, 5);
	$total_conect = substr(($_), 61, 5);
	$total_seqres = substr(($_), 66, 5);
    }

    #=============================================================================
    #   Keyword MODEL

    if ($first_field eq 'MODEL') {
	$xlat_flag = $xlat_save;
	$model_flag = $Fld[2];
	$model_count++;
	$model_flags = 'yes';
	$model_asn_low{$model_count} = $atom_flag;
	$model_asn_high{$model_count} = $atom_flag;
	$model_an_low{$model_count} = 0;
	$model_an_high{$model_count} = 0;
	if ($model_flag eq '') {
	    $model_flag = '.';
	}
	$model{$model_count} = $model_flag;
    }

    #===========================================================================
    #  Keyword MODRES
    #
    # In the 1995 format, a new record, MODRES, was added to provide
    # "descriptions of modifications (e.g., chemical or post-
    # translational) to protein and nucleic acid residues.  Inlcuded
    # are mapping between residue names given in a PDB entry
    # and standard residues."  We treat this record as if it were
    # a SEQADV with no database specified.  To complete the necessary
    # mmCIf category relationships, a dummy DBREF is created for
    # each chain involved.
    #
    #     MODRES           [1- 6]
    #     modres_idcode    [8-11]       = idcode of this entry (not used)
    #     modres_resName  [13-15]       = _struct_ref_seq_dif.mon_id
    #     modres_chainID   [ 17 ]       = _struct_ref.biol_id
    #     modres_seq      [19-23]       = combines seqNum and insertCode
    #                                     used to derive
    #                                     _struct_ref_seq_dif.seq_num 
    #     "."                           = _struct_ref.db_name
    #     "."
    #                                   = _struct_ref.db_code
    #     modres_dbRes    [25-27]       = _struct_ref_seq_dif.db_mon_id
    #     "."                           = _struct_ref_seq_dif.db_seq_num
    #     modres_conflict [30-70]       = _struct_ref_seq_dif.details

    if ($first_field eq 'MODRES') {
	$xlat_flag = $xlat_save;
	$modres_resName{++$modres_flag} = substr(($_), 13, 3);
	$modres_chainID{$modres_flag} = substr(($_), 17, 1);
	if ($modres_chainID{$modres_flag} eq ' ') {
	    $modres_chainID{$modres_flag} = $bcid;
	}
	$modres_seq{$modres_flag} = substr(($_), 19, 5);
	$modres_dbRes{$modres_flag} = substr(($_), 25, 3);
	$modres_conflict{$modres_flag} = substr(($_), 30, 41);
	if ($modres_dbSeq{$modres_flag} eq '      ') {
	    $modres_dbSeq{$modres_flag} = '.';
	}
	$modres_conflict{$modres_flag} = ('Chain ' .

	  $modres_chainID{$modres_flag} . ': ' .

	  $modres_conflict{$modres_flag});
	{
	    $numx = (@dblist = split(' ', $modres_conflict{$modres_flag},

	      9999));
	    if ($numx >= 1 && $dblist[$numx] eq '' && ' ' eq ' ') {
		--$numx;
	    }
	}
	$modres_conflict{$modres_flag} = '.';
	if ($numx > 0) {
	    $modres_conflict{$modres_flag} = $dblist[1];
	    for ($j = 2; $j <= $numx; ++$j) {
		$modres_conflict{$modres_flag} =

		  ($modres_conflict{$modres_flag} . ' ' . $dblist[$j]);
	    }
	    if ($numx > 1) {
		$modres_conflict{$modres_flag} = ("'" .

		  $modres_conflict{$modres_flag} . "'");
	    }
	}
    }

    #============================================================================
    #   Keyword MTRIX
    #

    if ($first_field eq 'MTRIX1' || $first_field eq 'MTRIX2' ||

      $first_field eq 'MTRIX3') {
	$xlat_flag = $xlat_save;

	$mtrix_col1 = substr(($_), 11, 10);
	$mtrix_col2 = substr(($_), 21, 10);
	$mtrix_col3 = substr(($_), 31, 10);
	$mtrix_col4 = substr(($_), 46, 10);
	# print loop headers

	if ($mtrix_flag eq '0') {	#???
	    $mat_save{++$mtrix_flag} = "\n\n\n";
	    $mat_save{++$mtrix_flag} = "##############################\n";
	    $mat_save{++$mtrix_flag} = "#                            #\n";
	    $mat_save{++$mtrix_flag} = "# STRUCT_NCS_OPER            #\n";
	    $mat_save{++$mtrix_flag} = "#                            #\n";
	    $mat_save{++$mtrix_flag} = "##############################\n";

	    $mat_save{++$mtrix_flag} = "\n";

	    $mat_save{++$mtrix_flag} =

	      "# **** WARNING ****  Domain information needed \n";
	    $warning_list{++$warning_flag} =

	      "#=# STRUCT_NCS_OPER: Domain information needed\n";

	    $mat_save{++$mtrix_flag} = "\nloop_ \n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.id\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.code\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.matrix[1][1]\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.matrix[1][2]\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.matrix[1][3]\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.vector[1] \n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.matrix[2][1]\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.matrix[2][2]\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.matrix[2][3]\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.vector[2] \n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.matrix[3][1]\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.matrix[3][2]\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.matrix[3][3]\n";
	    $mat_save{++$mtrix_flag} = "_struct_ncs_oper.vector[3] \n";
	}

	$mtrix_id = substr(($_), 8, 3);
	$mtrix_given = substr(($_), 60, 1);
	$x_given = 'generate';
	if ($mtrix_given ne ' ') {
	    $x_given = 'given';
	}
	if ($first_field eq 'MTRIX1') {
	    $mat_save{++$mtrix_flag} = sprintf("%3s    %s\n", $mtrix_id,

	      $x_given);
	}
	$mat_save{++$mtrix_flag} = ($mtrix_col1 . ' ' . $mtrix_col2 . ' ' .

	  $mtrix_col3 . ' ' . $mtrix_col4 . "\n");
    }

    #===========================================================================
    #  Keyword OBSLTE:  see SPRSDE, below 
    #
    #============================================================================
    #  Keyword ORIGX
    #
    # _database_pdb_matrix.origx[1][1] .. [3][3]
    # _database_pdb_matrix.origx_vector[1]  .. _3

    if ($first_field eq 'ORIGX1' || $first_field eq 'ORIGX2' ||

      $first_field eq 'ORIGX3') {
	$xlat_flag = $xlat_save;

	$origx_col1 = substr(($_), 11, 10);
	$origx_col2 = substr(($_), 21, 10);
	$origx_col3 = substr(($_), 31, 10);
	$origx_col4 = substr(($_), 46, 10);

	# print loop headers

	if ($origx_flag eq '0') {	#???
	    $om_save{++$origx_flag} = "\n\n\n";

	    $om_save{++$origx_flag} = "\nloop_ \n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.entry_id\n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.origx[1][1]\n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.origx[1][2]\n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.origx[1][3]\n";
	    $om_save{++$origx_flag} =

	      "_database_pdb_matrix.origx_vector[1] \n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.origx[2][1]\n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.origx[2][2]\n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.origx[2][3]\n";
	    $om_save{++$origx_flag} =

	      "_database_pdb_matrix.origx_vector[2] \n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.origx[3][1]\n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.origx[3][2]\n";
	    $om_save{++$origx_flag} = "_database_pdb_matrix.origx[3][3]\n";
	    $om_save{++$origx_flag} =

	      "_database_pdb_matrix.origx_vector[3] \n";
	    $om_save{++$origx_flag} = ('  ' . $head_PDB_code . "\n\n");
	}

	$origx_id = substr(($_), 8, 3);
	$om_save{++$origx_flag} = ($origx_col1 . ' ' . $origx_col2 . ' ' .

	  $origx_col3 . ' ' . $origx_col4 . "\n");
    }

    #===========================================================================
    #   Keyword REMARK
    #
    #
    #  print all citations from JNRL and REMARK 1 records
    #
    #  First check if it is time to flush references
    #
    if ($flush_refs == 1) {
	$remark_number = substr(($_), 8, 3);
	if (($first_field eq 'REMARK' && $remark_number eq '  2' &&	#???

	  $Fld[3] eq 'RESOLUTION.') ||

	  (($first_field ne 'REMARK') && ($first_field ne 'JRNL'))) {
	    if ($jrnl_flag eq '1') {	#???
		printf (("\nloop_\n"));
		printf (("_citation.id\n"));
		printf (("_citation.coordinate_linkage\n"));
		printf (("_citation.title\n"));
		printf (("_citation.country\n"));
		printf (("_citation.journal_abbrev\n"));
		printf (("_citation.journal_volume\n"));
		printf (("_citation.journal_issue\n"));
		printf (("_citation.page_first\n"));
		printf (("_citation.year\n"));
		printf (("_citation.journal_id_ASTM\n"));
		printf (("_citation.journal_id_ISSN\n"));
		printf (("_citation.journal_id_CSD\n"));
		printf (("_citation.book_title\n"));
		printf (("_citation.book_publisher\n"));
		printf (("_citation.book_id_ISBN\n"));
		printf (("_citation.details\n"));
		++$jrnl_flag;
	    }
	    $cit_decr = 0;
	    if ($primary) {
		$cit_decr = 1;
	    }
	    for ($i = 1; $i <= $cit_flag; ++$i) {
		if ($i eq '1' && $primary) {	#???
		    printf ((" \nprimary   yes\n"));
		}
		else {
		    printf " \n%3s       no\n", $i - $cit_decr;
		    if ($i - $cit_decr != $cit_refNum{$i}) {
			$warning_list{++$warning_flag} =

			  sprintf("#=# CITATION: Mismatch PDB refNum %s to id %s\n",

			  $cit_refNum{$i}, $i - $cit_decr);
		    }
		}

		# for books

		#       _citation.title	                   == TITL (if present)
		#       _citation.country		   == country[i]
		#       _citation.journal_abbrev           == ?
		#       _citation.journal_volume           == volu[i]
		#       _citation.journal_issue            == ?
		#       _citation.page_first               == ?
		#       _citation.year                     == year[i]
		#       _citation.journal_id_ASTM       == ?
		#       _citation.journal_id_ISSN       == ?
		#       _citation.journal_id_PDB        == ?
		#       _citation.book_title               == REF (jour_x)
		#       _citation.book_publisher           == jrnl_pub_pub_x[i]
		#       _citation.book_id_ISBN          == isbn[i]
		#       _citation.details          == ?

		if ($jrnl_pub_pub_1{$i} ne '?' || $issn_isbn{$i} eq 'ISBN') {
		    if ($country{$i} eq '  ' || $country{$i} eq '') {
			$country{$i} = '?';
		    }
		    if ($jour_1{$i} eq '                            ') {
			$jour_1{$i} = '?';
		    }
		    if ($volu{$i} eq '    ' || $volu{$i} eq '') {
			$volu{$i} = '?';
		    }
		    if ($year{$i} eq '    ' || $year{$i} eq '') {
			$year{$i} = '?';
		    }
		    if ($page{$i} eq '     ' || $page{$i} eq '') {
			$page{$i} = '?';
		    }
		    if ($csd{$i} eq '    ' || $csd{$i} eq '') {
			$csd{$i} = '?';
		    }
		    if ($cit_title_1{$i}) {
			printf "; %s\n", $cit_title_1{$i};
		    }
		    else {
			printf ((' ? '));
		    }
		    if ($cit_title_2{$i}) {
			printf "  %s\n", $cit_title_2{$i};
		    }
		    if ($cit_title_1{$i}) {
			printf ((";\n"));
		    }
		    printf " %2s ? %3s ? %5s %4s ? ? %4s\n", $country{$i},

		      $volu{$i}, $page{$i}, $year{$i}, $csd{$i};

		    if (!$jour_2{$i}) {
			printf " '%28s' \n", $jour_1{$i};
		    }

		    if ($jour_2{$i}) {
			printf "; %28s\n  %s\n;\n", $jour_1{$i}, $jour_2{$i};
		    }

		    if ($jrnl_pub_pub_1{$i}) {
			printf ";  %s \n", $jrnl_pub_pub_1{$i};
		    }
		    if ($jrnl_pub_pub_2{$i}) {
			printf "   %s \n", $jrnl_pub_pub_2{$i};
		    }
		    if ($jrnl_pub_pub_1{$i}) {
			printf ((";\n"));
		    }
		    printf " '%25s' ? \n", $isbn{$i};
		}

		else {
		    # for journals

		    if ($cit_title_1{$i}) {
			printf "; %s\n", $cit_title_1{$i};
		    }
		    if ($cit_title_2{$i}) {
			printf "  %s\n", $cit_title_2{$i};
		    }
		    if ($country{$i} eq '  ' || $country{$i} eq '') {
			$country{$i} = '?';
		    }
		    if ($volu{$i} eq '    ' || $volu{$i} eq '') {
			$volu{$i} = '?';
		    }
		    if ($year{$i} eq '    ' || $year{$i} eq '') {
			$year{$i} = '?';
		    }
		    if ($page{$i} eq '     ' || $page{$i} eq '') {
			$page{$i} = '?';
		    }
		    if ($csd{$i} eq '    ' || $csd{$i} eq '') {
			$csd{$i} = '?';
		    }
		    if (!$jour_2{$i}) {
			printf

			  ";\n %2s '%28s' %4s  ?   %5s   %4s \n'%-15s' '%15s' %4s ? ? ? ?\n",

			  $country{$i}, $jour_1{$i}, $volu{$i}, $page{$i},

			  $year{$i}, $astm{$i}, $issn{$i}, $csd{$i};
		    }
		    if ($jour_2{$i}) {
			printf ((";\n %2s \n; %-28s\n  %s\n;\n" .

			  " %4s  ?   %5s   %4s \n'%-15s' '%15s' %4s ? ? ? ?\n"),

			  $country{$i}, $jour_1{$i}, $jour_2{$i}, $volu{$i},

			  $page{$i}, $year{$i}, $astm{$i}, $issn{$i},

			  $csd{$i});
		    }
		}
	    }
	    # Loop Editor List

	    for ($i = 1; $i <= $cit_flag; ++$i) {
		if ($cit_edit_1{$i}) {
		    printf (("\nloop_\n"));
		    printf (("_citation_editor.citation_id\n"));
		    printf (("_citation_editor.name\n"));
		    last;
		}
	    }
	    for ($i = 1; $i <= $cit_flag; ++$i) {
		if ($cit_edit_1{$i}) {
		    {
			$num_edit = (@editors = split($comma, $cit_edit_1{$i},

			  9999));
			if ($num_edit >= 1 && $editors[$num_edit] eq '' &&

			  $comma eq ' ') {
			    --$num_edit;
			}
		    }
		    for ($ii = 1; $ii <= $num_edit; ++$ii) {
			{
			    $num_e_split = (@e_split = split(' ',

			      $editors[$ii], 9999));
			    if ($num_e_split >= 1 &&

			      $e_split[$num_e_split] eq '' && ' ' eq ' ') {
				--$num_e_split;
			    }
			}
			$editors[$ii] = '';
			if ($num_e_split > 0) {
			    $editors[$ii] = $e_split[1];
			    for ($j = 2; $j <= $num_e_split; ++$j) {
				$editors[$ii] = ($editors[$ii] . ' ' .

				  $e_split[$j]);
			    }
			}
			if ($auth_convtext eq 'yes' ||

			  ($auth_convtext eq 'conditional' &&

			  $convtext eq 'yes')) {
			    {
				#
				# produce a CIF-style name from a PDB name
				#
				# begin by applying typesetting codes if any
				# but always treat "-" and "'"  as breaks for capitalization
				# in names
				#
				{
				    $lx_tl = length($editors[$ii]);
				    $tx_tl = $editors[$ii];
				    $lostr = '';
				    for ($ix_tl = 1; $ix_tl <= $lx_tl;

				    ++$ix_tl) {
					$cx_tl = substr($tx_tl, $ix_tl, 1);
					$cx_tl = substr(($lcaz . $cx_tl),

					  index(($UCAZ . $cx_tl), $cx_tl), 1);
					$lostr = ($lostr . $cx_tl);
				    }
				}

				$lstr = length($lostr);
				$mystr = '';
				$pchar = ' ';
				for ($qnsi = 1; $qnsi <= $lstr; ++$qnsi) {
				    $mychar = substr($lostr, $qnsi, 1);
				    if ($pchar eq ' ' || $pchar eq ',' ||

				      $pchar eq '.' || $pchar eq '-' ||

				      $pchar eq "'" || $pchar eq '(' ||

				      $pchar eq '*' || $pchar eq '/') {
					$mychar = substr(($UCAZ . $mychar),

					  index(($lcaz . $mychar), $mychar),

					  1);
				    }
				    if (($mychar ne '*' && $mychar ne "\$" &&

				      $mychar ne '/') ||

				      ($mychar eq $pchar)) {	#???
					$mystr = ($mystr . $mychar);
				    }
				    if ($pchar eq '/') {
					if ($mychar eq "\$" ||

					  $mychar eq '-') {
					    $pchar = $mychar;
					}
					# end if( mychar == "$" || mychar == "-" ) 

					;
				    }
				    else {
					$pchar = $mychar;
				    }
				    # end if( pchar == "/" )

				    ;
				}
				# end for( qnsi=1; qnsi <= lstr; ++qnsi)
				#
				# See if a specific replacement was given
				#
				{
				    $lx_tu = length($mystr);
				    $tx_tu = $mystr;
				    $name_temp = '';
				    for ($ix_tu = 1; $ix_tu <= $lx_tu;

				    ++$ix_tu) {
					$cx_tu = substr($tx_tu, $ix_tu, 1);
					$cx_tu = substr(($UCAZ . $cx_tu),

					  index(($lcaz . $cx_tu), $cx_tu), 1);
					$name_temp = ($name_temp . $cx_tu);
				    }
				}

				if ($rep_name{$name_temp} ne '') {
				    $mystr = $rep_name{$name_temp};
				    #
				    # See if there is a comma in place if so we are done
				    #
				    ;
				}
				if (index($mystr, $comma) != 0) {
				    $ret_val = $mystr;
				}
				else {
				    $nam_suf = '';
				    {
					$num_namp = (@x_namep = split(' ',

					  $mystr, 9999));
					if ($num_namp >= 1 &&

					  $x_namep[$num_namp] eq '' &&

					  ' ' eq ' ') {
					    --$num_namp;
					}
				    }
				    if ($num_namp > 1) {
					{
					    $lx_tu =

					      length($x_namep[$num_namp]);
					    $tx_tu = $x_namep[$num_namp];
					    $xtemp = '';
					    for ($ix_tu = 1; $ix_tu <= $lx_tu;

					    ++$ix_tu) {
						$cx_tu = substr($tx_tu,

						  $ix_tu, 1);
						$cx_tu = substr(($UCAZ .

						  $cx_tu), index(($lcaz .

						  $cx_tu), $cx_tu), 1);
						$xtemp = ($xtemp . $cx_tu);
					    }
					}

					if ($rep_suffix{$xtemp} ne '') {
					    if ($junior_on_last eq 'yes') {
						$x_namep[$num_namp] =

						  $rep_suffix{$xtemp};
					    }
					    else {
						$nam_suf = (' ' .

						  $rep_suffix{$xtemp});
						--$num_namp;
					    }
					}
					$mystr = $x_namep[1];
					for ($knamp = 2; $knamp <= $num_namp;

					++$knamp) {
					    $mystr = ($mystr . ' ' .

					      $x_namep[$knamp]);
					}
				    }
				    # end if (num_namp > 1)
				    $llname = length($mystr);
				    $cc = '';
				    for ($kc = $llname - 1; $kc > 1; --$kc) {
					$cp = $cc;
					$cc = substr($mystr, $kc, 1);
					if ($cc eq '.') {	#???
					    if ($cp ne ' ') {
						$mystr = (substr($mystr,

						  $kc + 1, $llname - $kc) .

						  $comma . ' ' .

						  substr($mystr, 1, $kc));
					    }
					    else {
						$mystr = (substr($mystr,

						  $kc + 2, $llname - $kc - 1)

						  . $comma . ' ' .

						  substr($mystr, 1, $kc));
					    }
					    # if (cp != " ")
					    $kc = 0;
					}
					# end if (cc == ".")

					;
				    }
				    # for (kc=llname-1; kc>1; --kc)
				    $mystr = ($mystr . $nam_suf);
				    $ret_val = $mystr;
				}
				# end if (index(mystr,comma) != 0 )

				$editors[$ii] = $ret_val;
			    }
			}
		    }
		    if ($cit_edit_2{$i}) {
			($num_edit = $num_edit - 1);
		    }
		    for ($j = 1; $j <= $num_edit; ++$j) {
			if (($primary) && $i eq '1') {	#???
			    printf " primary   '%s' \n", $editors[$j];
			}
			else {
			    printf " %3s   '%s' \n", $i - $cit_decr,

			      $editors[$j];
			}
		    }
		}

		if ($cit_edit_2{$i}) {
		    {
			$num_edit = (@editors = split($comma, $cit_edit_2{$i},

			  9999));
			if ($num_edit >= 1 && $editors[$num_edit] eq '' &&

			  $comma eq ' ') {
			    --$num_edit;
			}
		    }
		    for ($ii = 1; $ii <= $num_edit; ++$ii) {
			{
			    $num_e_split = (@e_split = split(' ',

			      $editors[$ii], 9999));
			    if ($num_e_split >= 1 &&

			      $e_split[$num_e_split] eq '' && ' ' eq ' ') {
				--$num_e_split;
			    }
			}
			$editors[$ii] = '';
			if ($num_e_split > 0) {
			    $editors[$ii] = $e_split[1];
			    for ($j = 2; $j <= $num_e_split; ++$j) {
				$editors[$ii] = ($editors[$ii] . ' ' .

				  $e_split[$j]);
			    }
			}
			if ($auth_convtext eq 'yes' ||

			  ($auth_convtext eq 'conditional' &&

			  $convtext eq 'yes')) {
			    {
				#
				# produce a CIF-style name from a PDB name
				#
				# begin by applying typesetting codes if any
				# but always treat "-" and "'"  as breaks for capitalization
				# in names
				#
				{
				    $lx_tl = length($editors[$ii]);
				    $tx_tl = $editors[$ii];
				    $lostr = '';
				    for ($ix_tl = 1; $ix_tl <= $lx_tl;

				    ++$ix_tl) {
					$cx_tl = substr($tx_tl, $ix_tl, 1);
					$cx_tl = substr(($lcaz . $cx_tl),

					  index(($UCAZ . $cx_tl), $cx_tl), 1);
					$lostr = ($lostr . $cx_tl);
				    }
				}

				$lstr = length($lostr);
				$mystr = '';
				$pchar = ' ';
				for ($qnsi = 1; $qnsi <= $lstr; ++$qnsi) {
				    $mychar = substr($lostr, $qnsi, 1);
				    if ($pchar eq ' ' || $pchar eq ',' ||

				      $pchar eq '.' || $pchar eq '-' ||

				      $pchar eq "'" || $pchar eq '(' ||

				      $pchar eq '*' || $pchar eq '/') {
					$mychar = substr(($UCAZ . $mychar),

					  index(($lcaz . $mychar), $mychar),

					  1);
				    }
				    if (($mychar ne '*' && $mychar ne "\$" &&

				      $mychar ne '/') ||

				      ($mychar eq $pchar)) {	#???
					$mystr = ($mystr . $mychar);
				    }
				    if ($pchar eq '/') {
					if ($mychar eq "\$" ||

					  $mychar eq '-') {
					    $pchar = $mychar;
					}
					# end if( mychar == "$" || mychar == "-" ) 

					;
				    }
				    else {
					$pchar = $mychar;
				    }
				    # end if( pchar == "/" )

				    ;
				}
				# end for( qnsi=1; qnsi <= lstr; ++qnsi)
				#
				# See if a specific replacement was given
				#
				{
				    $lx_tu = length($mystr);
				    $tx_tu = $mystr;
				    $name_temp = '';
				    for ($ix_tu = 1; $ix_tu <= $lx_tu;

				    ++$ix_tu) {
					$cx_tu = substr($tx_tu, $ix_tu, 1);
					$cx_tu = substr(($UCAZ . $cx_tu),

					  index(($lcaz . $cx_tu), $cx_tu), 1);
					$name_temp = ($name_temp . $cx_tu);
				    }
				}

				if ($rep_name{$name_temp} ne '') {
				    $mystr = $rep_name{$name_temp};
				    #
				    # See if there is a comma in place if so we are done
				    #
				    ;
				}
				if (index($mystr, $comma) != 0) {
				    $ret_val = $mystr;
				}
				else {
				    $nam_suf = '';
				    {
					$num_namp = (@x_namep = split(' ',

					  $mystr, 9999));
					if ($num_namp >= 1 &&

					  $x_namep[$num_namp] eq '' &&

					  ' ' eq ' ') {
					    --$num_namp;
					}
				    }
				    if ($num_namp > 1) {
					{
					    $lx_tu =

					      length($x_namep[$num_namp]);
					    $tx_tu = $x_namep[$num_namp];
					    $xtemp = '';
					    for ($ix_tu = 1; $ix_tu <= $lx_tu;

					    ++$ix_tu) {
						$cx_tu = substr($tx_tu,

						  $ix_tu, 1);
						$cx_tu = substr(($UCAZ .

						  $cx_tu), index(($lcaz .

						  $cx_tu), $cx_tu), 1);
						$xtemp = ($xtemp . $cx_tu);
					    }
					}

					if ($rep_suffix{$xtemp} ne '') {
					    if ($junior_on_last eq 'yes') {
						$x_namep[$num_namp] =

						  $rep_suffix{$xtemp};
					    }
					    else {
						$nam_suf = (' ' .

						  $rep_suffix{$xtemp});
						--$num_namp;
					    }
					}
					$mystr = $x_namep[1];
					for ($knamp = 2; $knamp <= $num_namp;

					++$knamp) {
					    $mystr = ($mystr . ' ' .

					      $x_namep[$knamp]);
					}
				    }
				    # end if (num_namp > 1)
				    $llname = length($mystr);
				    $cc = '';
				    for ($kc = $llname - 1; $kc > 1; --$kc) {
					$cp = $cc;
					$cc = substr($mystr, $kc, 1);
					if ($cc eq '.') {	#???
					    if ($cp ne ' ') {
						$mystr = (substr($mystr,

						  $kc + 1, $llname - $kc) .

						  $comma . ' ' .

						  substr($mystr, 1, $kc));
					    }
					    else {
						$mystr = (substr($mystr,

						  $kc + 2, $llname - $kc - 1)

						  . $comma . ' ' .

						  substr($mystr, 1, $kc));
					    }
					    # if (cp != " ")
					    $kc = 0;
					}
					# end if (cc == ".")

					;
				    }
				    # for (kc=llname-1; kc>1; --kc)
				    $mystr = ($mystr . $nam_suf);
				    $ret_val = $mystr;
				}
				# end if (index(mystr,comma) != 0 )

				$editors[$ii] = $ret_val;
			    }
			}
		    }
		    for ($j = 1; $j <= $num_edit; ++$j) {
			if (($primary) && $i eq '1') {	#???
			    printf " primary   '%s' \n", $editors[$j];
			}
			else {
			    printf " %3s   '%s' \n", $i - $cit_decr,

			      $editors[$j];
			}
		    }
		}
	    }

	    # Loop Author List

	    for ($i = 1; $i <= $cit_flag; ++$i) {
		if ($cit_auth_1{$i}) {
		    printf (("\nloop_\n"));
		    printf (("_citation_author.citation_id\n"));
		    printf (("_citation_author.name\n"));
		    last;
		}
	    }

	    for ($i = 1; $i <= $cit_flag; ++$i) {
		if ($cit_auth_1{$i}) {
		    {
			$num_auth = (@authors = split($comma, $cit_auth_1{$i},

			  9999));
			if ($num_auth >= 1 && $authors[$num_auth] eq '' &&

			  $comma eq ' ') {
			    --$num_auth;
			}
		    }
		    for ($ii = 1; $ii <= $num_auth; ++$ii) {
			{
			    $num_a_split = (@a_split = split(' ',

			      $authors[$ii], 9999));
			    if ($num_a_split >= 1 &&

			      $a_split[$num_a_split] eq '' && ' ' eq ' ') {
				--$num_a_split;
			    }
			}
			$authors[$ii] = '';
			if ($num_a_split > 0) {
			    $authors[$ii] = $a_split[1];
			    for ($j = 2; $j <= $num_a_split; ++$j) {
				$authors[$ii] = ($authors[$ii] . ' ' .

				  $a_split[$j]);
			    }
			}
			if ($auth_convtext eq 'yes' ||

			  ($auth_convtext eq 'conditional' &&

			  $convtext eq 'yes')) {
			    {
				#
				# produce a CIF-style name from a PDB name
				#
				# begin by applying typesetting codes if any
				# but always treat "-" and "'"  as breaks for capitalization
				# in names
				#
				{
				    $lx_tl = length($authors[$ii]);
				    $tx_tl = $authors[$ii];
				    $lostr = '';
				    for ($ix_tl = 1; $ix_tl <= $lx_tl;

				    ++$ix_tl) {
					$cx_tl = substr($tx_tl, $ix_tl, 1);
					$cx_tl = substr(($lcaz . $cx_tl),

					  index(($UCAZ . $cx_tl), $cx_tl), 1);
					$lostr = ($lostr . $cx_tl);
				    }
				}

				$lstr = length($lostr);
				$mystr = '';
				$pchar = ' ';
				for ($qnsi = 1; $qnsi <= $lstr; ++$qnsi) {
				    $mychar = substr($lostr, $qnsi, 1);
				    if ($pchar eq ' ' || $pchar eq ',' ||

				      $pchar eq '.' || $pchar eq '-' ||

				      $pchar eq "'" || $pchar eq '(' ||

				      $pchar eq '*' || $pchar eq '/') {
					$mychar = substr(($UCAZ . $mychar),

					  index(($lcaz . $mychar), $mychar),

					  1);
				    }
				    if (($mychar ne '*' && $mychar ne "\$" &&

				      $mychar ne '/') ||

				      ($mychar eq $pchar)) {	#???
					$mystr = ($mystr . $mychar);
				    }
				    if ($pchar eq '/') {
					if ($mychar eq "\$" ||

					  $mychar eq '-') {
					    $pchar = $mychar;
					}
					# end if( mychar == "$" || mychar == "-" ) 

					;
				    }
				    else {
					$pchar = $mychar;
				    }
				    # end if( pchar == "/" )

				    ;
				}
				# end for( qnsi=1; qnsi <= lstr; ++qnsi)
				#
				# See if a specific replacement was given
				#
				{
				    $lx_tu = length($mystr);
				    $tx_tu = $mystr;
				    $name_temp = '';
				    for ($ix_tu = 1; $ix_tu <= $lx_tu;

				    ++$ix_tu) {
					$cx_tu = substr($tx_tu, $ix_tu, 1);
					$cx_tu = substr(($UCAZ . $cx_tu),

					  index(($lcaz . $cx_tu), $cx_tu), 1);
					$name_temp = ($name_temp . $cx_tu);
				    }
				}

				if ($rep_name{$name_temp} ne '') {
				    $mystr = $rep_name{$name_temp};
				    #
				    # See if there is a comma in place if so we are done
				    #
				    ;
				}
				if (index($mystr, $comma) != 0) {
				    $ret_val = $mystr;
				}
				else {
				    $nam_suf = '';
				    {
					$num_namp = (@x_namep = split(' ',

					  $mystr, 9999));
					if ($num_namp >= 1 &&

					  $x_namep[$num_namp] eq '' &&

					  ' ' eq ' ') {
					    --$num_namp;
					}
				    }
				    if ($num_namp > 1) {
					{
					    $lx_tu =

					      length($x_namep[$num_namp]);
					    $tx_tu = $x_namep[$num_namp];
					    $xtemp = '';
					    for ($ix_tu = 1; $ix_tu <= $lx_tu;

					    ++$ix_tu) {
						$cx_tu = substr($tx_tu,

						  $ix_tu, 1);
						$cx_tu = substr(($UCAZ .

						  $cx_tu), index(($lcaz .

						  $cx_tu), $cx_tu), 1);
						$xtemp = ($xtemp . $cx_tu);
					    }
					}

					if ($rep_suffix{$xtemp} ne '') {
					    if ($junior_on_last eq 'yes') {
						$x_namep[$num_namp] =

						  $rep_suffix{$xtemp};
					    }
					    else {
						$nam_suf = (' ' .

						  $rep_suffix{$xtemp});
						--$num_namp;
					    }
					}
					$mystr = $x_namep[1];
					for ($knamp = 2; $knamp <= $num_namp;

					++$knamp) {
					    $mystr = ($mystr . ' ' .

					      $x_namep[$knamp]);
					}
				    }
				    # end if (num_namp > 1)
				    $llname = length($mystr);
				    $cc = '';
				    for ($kc = $llname - 1; $kc > 1; --$kc) {
					$cp = $cc;
					$cc = substr($mystr, $kc, 1);
					if ($cc eq '.') {	#???
					    if ($cp ne ' ') {
						$mystr = (substr($mystr,

						  $kc + 1, $llname - $kc) .

						  $comma . ' ' .

						  substr($mystr, 1, $kc));
					    }
					    else {
						$mystr = (substr($mystr,

						  $kc + 2, $llname - $kc - 1)

						  . $comma . ' ' .

						  substr($mystr, 1, $kc));
					    }
					    # if (cp != " ")
					    $kc = 0;
					}
					# end if (cc == ".")

					;
				    }
				    # for (kc=llname-1; kc>1; --kc)
				    $mystr = ($mystr . $nam_suf);
				    $ret_val = $mystr;
				}
				# end if (index(mystr,comma) != 0 )

				$authors[$ii] = $ret_val;
			    }
			}
		    }
		    if ($cit_auth_2{$i}) {
			($num_auth = $num_auth - 1);
		    }
		    for ($j = 1; $j <= $num_auth; ++$j) {
			if (($primary) && $i eq '1') {	#???
			    printf " primary   '%s' \n", $authors[$j];
			}
			else {
			    printf " %3s       '%s' \n", $i - $cit_decr,

			      $authors[$j];
			}
		    }
		}

		if ($cit_auth_2{$i}) {
		    {
			$num_auth = (@authors = split($comma, $cit_auth_2{$i},

			  9999));
			if ($num_auth >= 1 && $authors[$num_auth] eq '' &&

			  $comma eq ' ') {
			    --$num_auth;
			}
		    }
		    for ($ii = 1; $ii <= $num_auth; ++$ii) {
			{
			    $num_a_split = (@a_split = split(' ',

			      $authors[$ii], 9999));
			    if ($num_a_split >= 1 &&

			      $a_split[$num_a_split] eq '' && ' ' eq ' ') {
				--$num_a_split;
			    }
			}
			$authors[$ii] = '';
			if ($num_a_split > 0) {
			    $authors[$ii] = $a_split[1];
			    for ($j = 2; $j <= $num_a_split; ++$j) {
				$authors[$ii] = ($authors[$ii] . ' ' .

				  $a_split[$j]);
			    }
			}
			if ($auth_convtext eq 'yes' ||

			  ($auth_convtext eq 'conditional' &&

			  $convtext eq 'yes')) {
			    {
				#
				# produce a CIF-style name from a PDB name
				#
				# begin by applying typesetting codes if any
				# but always treat "-" and "'"  as breaks for capitalization
				# in names
				#
				{
				    $lx_tl = length($authors[$ii]);
				    $tx_tl = $authors[$ii];
				    $lostr = '';
				    for ($ix_tl = 1; $ix_tl <= $lx_tl;

				    ++$ix_tl) {
					$cx_tl = substr($tx_tl, $ix_tl, 1);
					$cx_tl = substr(($lcaz . $cx_tl),

					  index(($UCAZ . $cx_tl), $cx_tl), 1);
					$lostr = ($lostr . $cx_tl);
				    }
				}

				$lstr = length($lostr);
				$mystr = '';
				$pchar = ' ';
				for ($qnsi = 1; $qnsi <= $lstr; ++$qnsi) {
				    $mychar = substr($lostr, $qnsi, 1);
				    if ($pchar eq ' ' || $pchar eq ',' ||

				      $pchar eq '.' || $pchar eq '-' ||

				      $pchar eq "'" || $pchar eq '(' ||

				      $pchar eq '*' || $pchar eq '/') {
					$mychar = substr(($UCAZ . $mychar),

					  index(($lcaz . $mychar), $mychar),

					  1);
				    }
				    if (($mychar ne '*' && $mychar ne "\$" &&

				      $mychar ne '/') ||

				      ($mychar eq $pchar)) {	#???
					$mystr = ($mystr . $mychar);
				    }
				    if ($pchar eq '/') {
					if ($mychar eq "\$" ||

					  $mychar eq '-') {
					    $pchar = $mychar;
					}
					# end if( mychar == "$" || mychar == "-" ) 

					;
				    }
				    else {
					$pchar = $mychar;
				    }
				    # end if( pchar == "/" )

				    ;
				}
				# end for( qnsi=1; qnsi <= lstr; ++qnsi)
				#
				# See if a specific replacement was given
				#
				{
				    $lx_tu = length($mystr);
				    $tx_tu = $mystr;
				    $name_temp = '';
				    for ($ix_tu = 1; $ix_tu <= $lx_tu;

				    ++$ix_tu) {
					$cx_tu = substr($tx_tu, $ix_tu, 1);
					$cx_tu = substr(($UCAZ . $cx_tu),

					  index(($lcaz . $cx_tu), $cx_tu), 1);
					$name_temp = ($name_temp . $cx_tu);
				    }
				}

				if ($rep_name{$name_temp} ne '') {
				    $mystr = $rep_name{$name_temp};
				    #
				    # See if there is a comma in place if so we are done
				    #
				    ;
				}
				if (index($mystr, $comma) != 0) {
				    $ret_val = $mystr;
				}
				else {
				    $nam_suf = '';
				    {
					$num_namp = (@x_namep = split(' ',

					  $mystr, 9999));
					if ($num_namp >= 1 &&

					  $x_namep[$num_namp] eq '' &&

					  ' ' eq ' ') {
					    --$num_namp;
					}
				    }
				    if ($num_namp > 1) {
					{
					    $lx_tu =

					      length($x_namep[$num_namp]);
					    $tx_tu = $x_namep[$num_namp];
					    $xtemp = '';
					    for ($ix_tu = 1; $ix_tu <= $lx_tu;

					    ++$ix_tu) {
						$cx_tu = substr($tx_tu,

						  $ix_tu, 1);
						$cx_tu = substr(($UCAZ .

						  $cx_tu), index(($lcaz .

						  $cx_tu), $cx_tu), 1);
						$xtemp = ($xtemp . $cx_tu);
					    }
					}

					if ($rep_suffix{$xtemp} ne '') {
					    if ($junior_on_last eq 'yes') {
						$x_namep[$num_namp] =

						  $rep_suffix{$xtemp};
					    }
					    else {
						$nam_suf = (' ' .

						  $rep_suffix{$xtemp});
						--$num_namp;
					    }
					}
					$mystr = $x_namep[1];
					for ($knamp = 2; $knamp <= $num_namp;

					++$knamp) {
					    $mystr = ($mystr . ' ' .

					      $x_namep[$knamp]);
					}
				    }
				    # end if (num_namp > 1)
				    $llname = length($mystr);
				    $cc = '';
				    for ($kc = $llname - 1; $kc > 1; --$kc) {
					$cp = $cc;
					$cc = substr($mystr, $kc, 1);
					if ($cc eq '.') {	#???
					    if ($cp ne ' ') {
						$mystr = (substr($mystr,

						  $kc + 1, $llname - $kc) .

						  $comma . ' ' .

						  substr($mystr, 1, $kc));
					    }
					    else {
						$mystr = (substr($mystr,

						  $kc + 2, $llname - $kc - 1)

						  . $comma . ' ' .

						  substr($mystr, 1, $kc));
					    }
					    # if (cp != " ")
					    $kc = 0;
					}
					# end if (cc == ".")

					;
				    }
				    # for (kc=llname-1; kc>1; --kc)
				    $mystr = ($mystr . $nam_suf);
				    $ret_val = $mystr;
				}
				# end if (index(mystr,comma) != 0 )

				$authors[$ii] = $ret_val;
			    }
			}
		    }
		    for ($j = 1; $j <= $num_auth; ++$j) {
			if (($primary) && $i eq '1') {	#???
			    printf " primary   '%s' \n", $authors[$j];
			}
			else {
			    printf " %3s       '%s' \n", $i - $cit_decr,

			      $authors[$j];
			}
		    }
		}
	    }
	    $flush_refs = 0;
	}
    }

    if ($first_field eq 'REMARK') {
	$xlat_flag = $xlat_save;

	++$all_remarks;

	# parse record

	$remark_number = substr(($_), 8, 3);
	$remark_cont = substr(($_), 17, 2);
	$jrnl_rec_type = substr(($_), 13, 4);
	$jrnl_refNum = substr(($_), 22, 49) + 0;
	$remark_text = substr(($_), 12, 60);
	$remark_cit_text = substr(($_), 20, 51);
	if ($convtext eq 'yes' && $jrnl_rec_type eq 'TITL') {
	    {
		#
		# apply PDB typsetting codes if any to a line
		#
		{
		    $lx_tl = length($remark_cit_text);
		    $tx_tl = $remark_cit_text;
		    $lostr = '';
		    for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
			$cx_tl = substr($tx_tl, $ix_tl, 1);
			$cx_tl = substr(($lcaz . $cx_tl),

			  index(($UCAZ . $cx_tl), $cx_tl), 1);
			$lostr = ($lostr . $cx_tl);
		    }
		}

		$lstr = length($lostr);
		$mystr = '';
		$pchar = ' ';
		for ($qtsi = 1; $qtsi <= $lstr; ++$qtsi) {
		    $mychar = substr($lostr, $qtsi, 1);
		    if ($pchar eq ' ' || $pchar eq ',' || $pchar eq '.' ||

		      $pchar eq '(' || $pchar eq '*' || $pchar eq '/') {
			$mychar = substr(($UCAZ . $mychar),

			  index(($lcaz . $mychar), $mychar), 1);
		    }
		    if (($mychar ne '*' && $mychar ne "\$" &&

		      $mychar ne '/') || ($mychar eq $pchar)) {	#???
			$mystr = ($mystr . $mychar);
		    }
		    if ($pchar eq '/') {
			if ($mychar eq "\$" || $mychar eq '-') {
			    $pchar = $mychar;
			}
		    }
		    else {
			$pchar = $mychar;
		    }
		}
		$ret_val = $mystr;

		$remark_cit_text = $ret_val;
	    }
	}
	$remark_test = substr(($_), 12, 3);

	# Deal with change of remark number

	$remark_test = substr(($_), 12, 3);

	if ($remark_number_old != $remark_number) {	#???
	    $remark_flag = '0';
	    $remark_number_old = $remark_number;
	    if ($remark_number ne '  3' && $remark_number ne '  2' &&	#???	#???

	      $remark_number ne '  1') {	#???
		printf ((";\n\n"));
	    }
	}

	#
	# As of the February 1996 PDB format, Remark 4 contains text
	# indicating the format with which the entry complies
	#
	if ($remark_number eq '  4') {	#???
	    if (substr($_, 17, 23) eq 'COMPLIES WITH FORMAT V.') {
		$compliance_level = substr($_, 41, 3);
	    }
	}
	# type 1 remarks - additional references
	# data items identical to JRNL

	if ($remark_number eq '  1' && $remark_test eq 'REF') {	#???
	    ++$cit_flag;
	    $flush_refs = 1;
	    $cit_refNum{$cit_flag} = $jrnl_refNum;
	}
	if ($remark_number eq '  1' && $remark_test ne 'REF' &&	#???

	  $remark_test ne '   ') {
	    # Assign TITL records

	    if ($jrnl_rec_type eq 'TITL' && $remark_cont eq '  ') {
		$cit_title_1{$cit_flag} = $remark_cit_text;
		$cit_title_2{$cit_flag} = '';
	    }
	    if ($jrnl_rec_type eq 'TITL' && $remark_cont ne '  ') {
		if ($cit_title_2{$cit_flag} eq '') {
		    $cit_title_2{$cit_flag} = $remark_cit_text;
		}
		else {
		    $cit_title_2{$cit_flag} = ($cit_title_2{$cit_flag} .

		      "\n  " . $remark_cit_text);
		}
	    }

	    # Assign AUTH records

	    if ($jrnl_rec_type eq 'AUTH' && $remark_cont eq '  ') {
		$cit_auth_1{$cit_flag} = $remark_cit_text;
		$cit_auth_2{$cit_flag} = '';
	    }
	    if ($jrnl_rec_type eq 'AUTH' && $remark_cont ne '  ') {
		$cit_auth_2{$cit_flag} = ($cit_auth_2{$cit_flag} .

		  $remark_cit_text);

		# Assign EDIT records

		;
	    }
	    if ($jrnl_rec_type eq 'EDIT' && $remark_cont eq '  ') {
		$cit_edit_1{$cit_flag} = $remark_cit_text;
		$cit_edit_2{$cit_flag} = '';
	    }
	    if ($jrnl_rec_type eq 'EDIT' && $remark_cont ne '  ') {
		$cit_edit_2{$cit_flag} = ($cit_edit_2{$cit_flag} .

		  $remark_cit_text);

		# Assign REF records

		;
	    }
	    if ($jrnl_rec_type eq 'REF ' && $remark_cont eq '  ') {
		$jour_1{$cit_flag} = substr(($_), 20, 28);
		$jour_2{$cit_flag} = '';
		$volu{$cit_flag} = substr(($_), 52, 4);
		$page{$cit_flag} = substr(($_), 57, 5);
		$year{$cit_flag} = substr(($_), 63, 4);
		$jrnl_pub_pub_1{$cit_flag} = '?';
	    }
	    if ($jrnl_rec_type eq 'REF ' && $remark_cont ne '  ') {
		if ($jour_2{$cit_flag} eq '') {
		    $jour_2{$cit_flag} = substr(($_), 20, 28);
		}
		else {
		    $jour_2{$cit_flag} = ($jour_2{$cit_flag} . "\n  " .

		      substr(($_), 20, 28));
		}
	    }
	    # Assign PUBL records

	    if ($jrnl_rec_type eq 'PUBL' && $remark_cont eq '  ') {
		$jrnl_pub_pub_1{$cit_flag} = substr(($_), 20, 51);
		$jrnl_pub_pub_2{$cit_flag} = '';
	    }
	    if ($jrnl_rec_type eq 'PUBL' && $remark_cont ne '  ') {
		if ($jrnl_pub_pub_2{$cit_flag} eq '') {
		    $jrnl_pub_pub_2{$cit_flag} = substr(($_), 20, 51);
		}
		else {
		    $jrnl_pub_pub_2{$cit_flag} = ($jrnl_pub_pub_2{$cit_flag} .

		      "\n  " . substr(($_), 20, 51));
		}
	    }
	    if ($jrnl_rec_type eq 'REFN') {
		$astm{$cit_flag} = substr(($_), 25, 6);
		$country{$cit_flag} = substr(($_), 33, 2);
		if ($country{$cit_flag} eq '  ') {
		    $country{$cit_flag} = '?';
		}
		$issn_isbn{$cit_flag} = substr(($_), 36, 4);
		if ($issn_isbn{$cit_flag} eq '    ' &&

		  substr($jour_1{$cit_flag}, 1, 9) ne 'TO BE PUB') {
		    if ($jrnl_pub_pub_1{$cit_flag} ne '?') {
			$issn_isbn{$cit_flag} = 'ISBN';
		    }
		    if ($volu{$cit_flag} eq '' || $volu{$cit_flag} eq '?' ||

		      $volu{$cit_flag} eq '    ') {
			$issn_isbn{$cit_flag} = 'ISBN';
		    }
		}
		if ($issn_isbn{$cit_flag} ne 'ISBN') {
		    $isbn{$cit_flag} = '?';
		    $issn{$cit_flag} = substr(($_), 41, 25);
		}
		else {
		    $issn{$cit_flag} = '?';
		    $isbn{$cit_flag} = substr(($_), 41, 25);
		}
		$csd{$cit_flag} = substr(($_), 67, 4);
		if ($csd{$cit_flag} eq '    ') {
		    $csd{$cit_flag} = '?';
		}
	    }
	    ++$remark_flag;
	}
	#
	#   type 2 remarks - resolution
	#

	if ($remark_number eq '  2' && $Fld[3] eq 'RESOLUTION.') {	#???
	    $resolution = substr(($_), 23, 45);

	    {
		$num_split = (@res_split = split(' ', $resolution, 9999));
		if ($num_split >= 1 && $res_split[$num_split] eq '' &&

		  ' ' eq ' ') {
		    --$num_split;
		}
	    }
	    if ($res_split[1] ne 'NOT') {
		$res_flag = 1;
		printf "\n_reflns.entry_id             %s \n", $head_PDB_code;
		printf "_reflns.d_resolution_high    %8.2g \n", $res_split[1];
	    }
	    ++$remark_flag;

	    #   Include _exptl templates

	    if ($verbose eq 'yes') {
		printf (("_exptl.absorpt_coefficient_mu         ? \n"));
		printf (("_exptl.absorpt_correction_T_max       ? \n"));
		printf (("_exptl.absorpt_correction_type        ? \n"));
		printf (("_exptl.absorpt_process_details        ? \n\n"));

		printf (("_exptl_crystal.colour                 ? \n"));
		printf (("_exptl_crystal.density_diffrn         ? \n"));
		printf (("_exptl_crystal.density_meas           ? \n"));
		printf (("_exptl_crystal.density_meas_temp      ? \n"));
		printf (("_exptl_crystal.density_method         ? \n"));
		printf (("_exptl_crystal.description            ? \n"));
		printf (("_exptl_crystal.F_000                  ? \n"));
		printf (("_exptl_crystal_face.diffr_chi         ? \n"));
		printf (("_exptl_crystal_face.diffr_kappa       ? \n"));
		printf (("_exptl_crystal_face.diffr_phi         ? \n"));
		printf (("_exptl_crystal_face.diffr_psi         ? \n"));
		printf (("_exptl_crystal_face.index_h           ? \n"));
		printf (("_exptl_crystal_face.index_k           ? \n"));
		printf (("_exptl_crystal_face.index_l           ? \n"));
		printf (("_exptl_crystal_face.perp_dist         ? \n"));
		printf (("_exptl_crystal.id                     ? \n"));
		printf (("_exptl_crystal.preparation            ? \n"));
		printf (("_exptl_crystal.size_max               ? \n"));
		printf (("_exptl_crystal.size_mid               ? \n"));
		printf (("_exptl_crystal.size_min               ? \n"));
		printf (("_exptl_crystal.size_rad               ? \n"));
		printf (("_exptl.crystals_number                ? \n"));

		printf (("_exptl_crystal_grow.apparatus         ? \n"));
		printf (("_exptl_crystal_grow.atmosphere        ? \n"));
		printf (("_exptl_crystal_grow.crystal_id        ? \n"));
		printf (("_exptl_crystal_grow.details           ? \n"));
		printf (("_exptl_crystal_grow.method            ? \n"));
		printf (("_exptl_crystal_grow.method_ref        ? \n"));
		printf (("_exptl_crystal_grow.pH                ? \n"));
		printf (("_exptl_crystal_grow.pressure          ? \n"));
		printf (("_exptl_crystal_grow.seeding           ? \n"));
		printf (("_exptl_crystal_grow.seeding_ref       ? \n"));
		printf (("_exptl_crystal_grow.temp              ? \n"));
		printf (("_exptl_crystal_grow.time              ? \n"));

		printf (("\nloop_\n"));
		printf (("_exptl_crystal_grow_comp.crystal_id \n"));
		printf (("_exptl_crystal_grow_comp.id \n"));
		printf (("_exptl_crystal_grow_comp.conc \n"));
		printf (("_exptl_crystal_grow_comp.details \n"));
		printf (("_exptl_crystal_grow_comp.name \n"));
		printf (("_exptl_crystal_grow_comp.sol_id \n"));
		printf (("_exptl_crystal_grow_comp.volume \n"));
		printf ((" ?  ?  ?  ?  ?  ?  ? \n\n"));
	    }

	    #   Include additional data items to be added on diffraction experiment.
	    #   A rigourous treatment of REMARK 3 might be able to parse some of
	    #   this info.
	    if ($verbose eq 'yes') {
		printf (("_diffrn.ambient_temp       ? \n"));
		printf (("_diffrn.ambient_pressure          ? \n"));
		printf (("_diffrn_attenuator.code           ? \n"));
		printf (("_diffrn_attenuator.scale          ? \n"));
		printf (("_diffrn.details                   ? \n\n"));

		printf (("_diffrn.ambient_environment       ? \n"));
		printf (("_diffrn.crystal_support           ? \n"));
		printf (("_diffrn.crystal_treatment         ? \n\n"));

		printf (("_diffrn_measurement.method        ? \n"));
		printf (("_diffrn_measurement.details       ? \n"));
		printf (("_diffrn_measurement.device        ? \n"));
		printf (("_diffrn_measurement.device_details    ? \n"));
		printf (("_diffrn_measurement.device_type       ? \n"));

		printf (("_diffrn_orient_matrix.type        ? \n"));
		printf (("_diffrn_orient_matrix.UB[1][1]       ? \n"));
		printf (("_diffrn_orient_matrix.UB[1][2]       ? \n"));
		printf (("_diffrn_orient_matrix.UB[1][3]       ? \n"));
		printf (("_diffrn_orient_matrix.UB[2][1]       ? \n"));
		printf (("_diffrn_orient_matrix.UB[2][2]       ? \n"));
		printf (("_diffrn_orient_matrix.UB[2][3]       ? \n"));
		printf (("_diffrn_orient_matrix.UB[3][1]       ? \n"));
		printf (("_diffrn_orient_matrix.UB[3][2]       ? \n"));
		printf (("_diffrn_orient_matrix.UB[3][3]       ? \n\n"));

		printf (("loop_\n"));
		printf (("_diffrn_orient_refln.index_h\n"));
		printf (("_diffrn_orient_refln.index_k\n"));
		printf (("_diffrn_orient_refln.index_l\n"));
		printf (("_diffrn_orient_refln.angle_chi\n"));
		printf (("_diffrn_orient_refln.angle_kappa\n"));
		printf (("_diffrn_orient_refln.angle_phi\n"));
		printf (("_diffrn_orient_refln.angle_psi\n"));
		printf ((" ?   ?   ?  ?  ?  ?  ?\n\n"));

		printf (("_diffrn_radiation.filter_edge     ? \n"));
		printf (("_diffrn_radiation.inhomogeneity   ? \n"));
		printf (("_diffrn_radiation.monochromator   ? \n"));
		printf (("_diffrn_radiation.polarisn_norm   ? \n"));
		printf (("_diffrn_radiation.polarisn_ratio  ? \n"));
		printf (("_diffrn_radiation.collimation     ? \n"));
		printf (("_diffrn_radiation.type            ? \n\n"));

		printf (("loop_\n"));
		printf (("_diffrn_radiation_wavelength.id \n"));
		printf (("_diffrn_radiation_wavelength.wavelength \n"));
		printf (("_diffrn_radiation_wavelength.wt \n"));
		printf ((" ?   ?   ?  \n\n"));

		printf (("_diffrn_detector.detector  ? \n"));
		printf (("_diffrn_detector.dtime     ? \n"));
		printf (("_diffrn_detector.details   ? \n"));
		printf (("_diffrn_detector.type      ? \n\n"));

		printf (("_diffrn_source.source      ? \n"));
		printf (("_diffrn_source.current     ? \n"));
		printf (("_diffrn_source.details     ? \n"));
		printf (("_diffrn_source.diffrn_id   ? \n"));
		printf (("_diffrn_source.power       ? \n"));
		printf (("_diffrn_source.size        ? \n"));
		printf (("_diffrn_source.target      ? \n"));
		printf (("_diffrn_source.type        ? \n\n"));
		printf (("_diffrn_source.voltage     ? \n\n"));

		printf (("loop_\n"));
		printf (("_diffrn_refln.index_h \n"));
		printf (("_diffrn_refln.index_k \n"));
		printf (("_diffrn_refln.index_l \n"));
		printf (("_diffrn_refln.angle_chi \n"));
		printf (("_diffrn_refln.angle_kappa \n"));
		printf (("_diffrn_refln.angle_omega \n"));
		printf (("_diffrn_refln.angle_phi \n"));
		printf (("_diffrn_refln.angle_psi \n"));
		printf (("_diffrn_refln.angle_theta \n"));
		printf (("_diffrn_refln.attenuator_code \n"));
		printf (("_diffrn_refln.counts_bg_1 \n"));
		printf (("_diffrn_refln.counts_bg_2 \n"));
		printf (("_diffrn_refln.counts_net \n"));
		printf (("_diffrn_refln.counts_peak \n"));
		printf (("_diffrn_refln.counts_total \n"));
		printf (("_diffrn_refln.detect_slit_horiz \n"));
		printf (("_diffrn_refln.detect_slit_vert \n"));
		printf (("_diffrn_refln.diffrn_id \n"));
		printf (("_diffrn_refln.elapsed_time \n"));
		printf (("_diffrn_refln.intensity_net \n"));
		printf (("_diffrn_refln.intensity_sigma \n"));
		printf (("_diffrn_refln.scale_group_code \n"));
		printf (("_diffrn_refln.scan_mode \n"));
		printf (("_diffrn_refln.scan_mode_backgd \n"));
		printf (("_diffrn_refln.scan_width \n"));
		printf (("_diffrn_refln.sint_over_lambda \n"));
		printf (("_diffrn_refln.standard_code\n"));
		printf (("_diffrn_refln.wavelength \n"));
		printf (("_diffrn_refln.wavelength_id \n"));
		printf

		  ((" ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ?\n\n"));

		printf (("_diffrn_reflns.av_R_equivalents      ? \n"));
		printf (("_diffrn_reflns.av_sigmaI_over_netI   ? \n"));
		printf (("_diffrn_reflns.limit_h_max           ? \n"));
		printf (("_diffrn_reflns.limit_h_min           ? \n"));
		printf (("_diffrn_reflns.limit_k_max           ? \n"));
		printf (("_diffrn_reflns.limit_k_min           ? \n"));
		printf (("_diffrn_reflns.limit_l_max           ? \n"));
		printf (("_diffrn_reflns.limit_l_min           ? \n"));
		printf (("_diffrn_reflns.number                ? \n"));
		printf (("_diffrn_reflns.reduction_process     ? \n"));
		printf (("_diffrn_reflns.theta_max             ? \n"));
		printf (("_diffrn_reflns.theta_min             ? \n"));
		printf (("_diffrn_reflns.transf_matrix[1][1]   ? \n"));
		printf (("_diffrn_reflns.transf_matrix[1][2]   ? \n"));
		printf (("_diffrn_reflns.transf_matrix[1][3]   ? \n"));
		printf (("_diffrn_reflns.transf_matrix[2][1]   ? \n"));
		printf (("_diffrn_reflns.transf_matrix[2][2]   ? \n"));
		printf (("_diffrn_reflns.transf_matrix[2][3]   ? \n"));
		printf (("_diffrn_reflns.transf_matrix[3][1]   ? \n"));
		printf (("_diffrn_reflns.transf_matrix[3][2]   ? \n"));
		printf (("_diffrn_reflns.transf_matrix[3][3]   ? \n\n"));

		printf (("loop_\n"));
		printf (("_diffrn_scale_group.code \n"));
		printf (("_diffrn_scale_group.I_net \n"));
		printf ((" ? ? \n\n"));

		printf (("loop_\n"));
		printf (("_diffrn_standard_refln.index_h \n"));
		printf (("_diffrn_standard_refln.index_k \n"));
		printf (("_diffrn_standard_refln.index_l \n"));
		printf (("_diffrn_standard_refln.code \n"));
		printf ((" ? ? ? ? \n\n"));

		printf (("_diffrn_standards.decay_%        ? \n"));
		printf (("_diffrn_standards.interval_count ? \n"));
		printf (("_diffrn_standards.interval_time  ? \n"));
		printf (("_diffrn_standards.number         ? \n"));
		printf (("_diffrn_standards.scale_sigma    ? \n"));

		printf (("\nloop_\n"));
		printf (("_refln.index_h \n"));
		printf (("_refln.index_k \n"));
		printf (("_refln.index_l \n"));
		printf (("_refln.A_meas \n"));
		printf (("_refln.A_calc \n"));
		printf (("_refln.B_meas \n"));
		printf (("_refln.B_calc \n"));
		printf (("_refln.crystal_id \n"));
		printf (("_refln.F_meas \n"));
		printf (("_refln.F_calc \n"));
		printf (("_refln.F_meas_sigma \n"));
		printf (("_refln.F_squared_meas \n"));
		printf (("_refln.F_squared_calc \n"));
		printf (("_refln.F_squared_sigma \n"));
		printf (("_refln.intensity_meas \n"));
		printf (("_refln.intensity_calc \n"));
		printf (("_refln.intensity_sigma \n"));
		printf (("_refln.mean_path_length_tbar \n"));
		printf (("_refln.status \n"));
		printf (("_refln.phase_meas \n"));
		printf (("_refln.phase_calc \n"));
		printf (("_refln.refinement_status \n"));
		printf (("_refln.scale_group_code \n"));
		printf (("_refln.sint_over_lambda \n"));
		printf (("_refln.symmetry_epsilon \n"));
		printf (("_refln.symmetry_multiplicity \n"));
		printf (("_refln.wavelength \n"));
		printf (("_refln.wavelength_id \n"));
		printf

		  ((" ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ?\n\n"));

		if ($res_flag == 0) {
		    printf "_reflns.entry_id              %s \n",

		      $head_PDB_code;
		    printf (("_reflns.d_resolution_high     ? \n"));
		}
		printf (("_reflns.d_resolution_low      ? \n"));
		printf (("_reflns.B_iso_Wilson_estimate ? \n"));
		printf (("_reflns.data_reduction_details ? \n"));
		printf (("_reflns.data_reduction_method ? \n"));
		printf (("_reflns.details              ? \n"));
		printf (("_reflns.limit_h_max           ? \n"));
		printf (("_reflns.limit_h_min           ? \n"));
		printf (("_reflns.limit_k_max           ? \n"));
		printf (("_reflns.limit_k_min           ? \n"));
		printf (("_reflns.limit_l_max           ? \n"));
		printf (("_reflns.limit_l_min           ? \n"));
		printf (("_reflns.number_all            ? \n"));
		printf (("_reflns.number_obs            ? \n"));
		printf (("_reflns.observed_criterion    ? \n"));
		printf (("_reflns.observed_criterion_sigma_F ? \n"));
		printf (("_reflns.observed_criterion_sigma_I ? \n"));
		printf (("_reflns.percent_possible_obs  ? \n"));
		printf (("_reflns.R_free_details        ? \n"));
		printf (("_reflns.Rmerge_F_all          ? \n"));
		printf (("_reflns.Rmerge_F_obs          ? \n"));

		printf (("\nloop_\n"));
		printf (("_reflns_scale.group_code \n"));
		printf (("_reflns_scale.meas_F \n"));
		printf (("_reflns_scale.meas_F_squared \n"));
		printf (("_reflns_scale.meas_intensity \n"));
		printf ((" ? ? ? ? \n\n"));

		printf (("_reflns.details               ? \n"));
		printf (("_reflns.data_reduction_method ? \n\n"));

		printf (("\nloop_\n"));
		printf (("_reflns_shell.d_res_high \n"));
		printf (("_reflns_shell.d_res_low  \n"));
		printf (("_reflns_shell.number_measured_all \n"));
		printf (("_reflns_shell.number_measured_obs \n"));
		printf (("_reflns_shell.number_possible \n"));
		printf (("_reflns_shell.number_unique_all \n"));
		printf (("_reflns_shell.number_unique_obs \n"));
		printf (("_reflns_shell.meanI_over_sigI_all \n"));
		printf (("_reflns_shell.meanI_over_sigI_obs \n"));
		printf (("_reflns_shell.percent_possible_all \n"));
		printf (("_reflns_shell.percent_possible_obs \n"));
		printf (("_reflns_shell.Rmerge_F_all \n"));
		printf (("_reflns_shell.Rmerge_I_all \n"));
		printf (("_reflns_shell.Rmerge_I_obs \n"));
		printf ((" ? ? ? ? ? ? ? ? ? ? ? ? ? ?\n\n"));

		printf (("_phasing_averaging.details       ? \n"));
		printf (("_phasing_averaging.method        ? \n"));

		printf (("_phasing_isomorphous.details     ? \n"));
		printf (("_phasing_isomorphous.method      ? \n"));
		printf (("_phasing_isomorphous.parent      ? \n"));

		printf (("_phasing_MAD.details             ? \n"));
		printf (("_phasing_MAD.method              ? \n"));

		printf (("_phasing_MIR.details             ? \n"));
		printf (("_phasing_MIR.method              ? \n"));

		printf (("\nloop_\n"));
		printf (("_phasing_MIR_der.id \n"));
		printf (("_phasing_MIR_der.number_of_sites \n"));
		printf (("_phasing_MIR_der.details \n"));
		printf (("_phasing_MIR_der.reflns_criteria \n"));
		printf ((" ?  ?  ?  ? \n"));

		printf (("\nloop_\n"));
		printf (("_phasing_MIR_der_shell.der_id \n"));
		printf (("_phasing_MIR_der_shell.d_res_low \n"));
		printf (("_phasing_MIR_der_shell.d_res_high \n"));
		printf (("_phasing_MIR_der_shell.fom \n"));
		printf (("_phasing_MIR_der_shell.ha_ampl \n"));
		printf (("_phasing_MIR_der_shell.loc \n"));
		printf (("_phasing_MIR_der_shell.phase \n"));
		printf (("_phasing_MIR_der_shell.power \n"));
		printf (("_phasing_MIR_der_shell.R_Cullis \n"));
		printf (("_phasing_MIR_der_shell.R_Kraut \n"));
		printf (("_phasing_MIR_der_shell.reflns \n"));
		printf ((" ?  ?  ?  ? ? ? ? ? ? ? ?\n"));

		printf (("\nloop_\n"));
		printf (("_phasing_MIR_shell.d_res_high \n"));
		printf (("_phasing_MIR_shell.d_res_low \n"));
		printf (("_phasing_MIR_shell.fom \n"));
		printf (("_phasing_MIR_shell.loc \n"));
		printf (("_phasing_MIR_shell.mean_phase \n"));
		printf (("_phasing_MIR_shell.power \n"));
		printf (("_phasing_MIR_shell.R_Cullis \n"));
		printf (("_phasing_MIR_shell.R_Kraut \n"));
		printf (("_phasing_MIR_shell.reflns \n"));
		printf ((" ?  ?  ?  ? ? ? ? ? ? \n"));

		printf (("\nloop_\n"));
		printf (("_phasing_MIR_der_site.der_id \n"));
		printf (("_phasing_MIR_der_site.id \n"));
		printf (("_phasing_MIR_der_site.B_iso \n"));
		printf (("_phasing_MIR_der_site.Cartn_x \n"));
		printf (("_phasing_MIR_der_site.Cartn_y \n"));
		printf (("_phasing_MIR_der_site.Cartn_z \n"));
		printf (("_phasing_MIR_der_site.fract_x \n"));
		printf (("_phasing_MIR_der_site.fract_y \n"));
		printf (("_phasing_MIR_der_site.fract_z \n"));
		printf (("_phasing_MIR_der_site.occupancy \n"));
		#       	printf ("_phasing_MIR_der_site.details \n")
		printf (("   ?  ?  ? ? ? ? ? ? ? ? \n\n"));
	    }
	}
	#
	#   type 3 remarks - refinmement details
	#

	if ($remark_number eq '  3') {	#???
	    if ($remark_flag eq '0') {	#???
		if ($verbose eq 'yes') {
		    printf (("_refine.diff_density_max          ? \n"));
		    printf (("_refine.diff_density_min          ? \n"));
		    printf (("_refine.ls_abs_structure_details  ? \n"));
		    printf (("_refine.ls_abs_structure_Flack    ? \n"));
		    printf (("_refine.ls_abs_structure_Rogers   ? \n"));
		    printf (("_refine.ls_extinction_coef        ? \n"));
		    printf (("_refine.ls_extinction_method      ? \n"));
		    printf (("_refine.ls_goodness_of_fit_all    ? \n"));
		    printf (("_refine.ls_goodness_of_fit_obs    ? \n"));
		    printf (("_refine.ls_hydrogen_treatment     ? \n"));
		    printf (("_refine.ls_matrix_type            ? \n"));
		    printf (("_refine.ls_number_constraints     ? \n"));
		    printf (("_refine.ls_number_parameters      ? \n"));
		    printf (("_refine.ls_number_reflns_obs      ? \n"));
		    printf (("_refine.ls_number_restraints      ? \n"));
		    printf (("_refine.ls_R_factor_all           ? \n"));
		    printf (("_refine.ls_R_factor_obs           ? \n"));
		    printf (("_refine.ls_restrained_S_all       ? \n"));
		    printf (("_refine.ls_restrained_S_obs       ? \n"));
		    printf (("_refine.ls_shift_over_esd_max     ? \n"));
		    printf (("_refine.ls_shift_over_esd_max     ? \n"));
		    printf (("_refine.ls_structure_factor_coef  ? \n"));
		    #       	printf ("_refine.ls_weighting_details      ? \n")
		    printf (("_refine.ls_weighting_scheme       ? \n"));
		    printf (("_refine.ls_wR_factor_all          ? \n"));
		    printf (("_refine.ls_wR_factor_obs          ? \n"));
		    printf (("_refine.details                   ? \n"));
		    printf (("_refine.occupancy_max             ? \n"));
		    printf (("_refine.occupancy_min             ? \n"));
		    printf (("_refine.B_iso_max                 ? \n"));
		    printf (("_refine.B_iso_min                 ? \n"));
		    printf (("_refine_ls_restr.criterion        ? \n"));
		    printf (("_refine_ls_restr.dev_ideal        ? \n"));
		    printf (("_refine_ls_restr.number           ? \n"));
		    printf (("_refine_ls_restr.rejects          ? \n"));
		    printf (("_refine_ls_restr.dev_ideal_target ? \n"));
		    printf (("_refine_ls_restr.type             ? \n"));

		    printf (("\nloop_\n"));
		    printf (("_refine_ls_shell.d_res_high \n"));
		    printf (("_refine_ls_shell.d_res_low  \n"));
		    printf (("_refine_ls_shell.number_reflns_all \n"));
		    printf (("_refine_ls_shell.number_reflns_obs \n"));
		    printf (("_refine_ls_shell.number_reflns_R_free \n"));
		    printf (("_refine_ls_shell.R_factor_all \n"));
		    printf (("_refine_ls_shell.R_factor_obs \n"));
		    printf (("_refine_ls_shell.wR_factor_all \n"));
		    printf (("_refine_ls_shell.wR_factor_obs \n"));
		    printf (("  ?  ?  ?  ?  ?  ?  ?  ?  ?\n"));

		    printf (("\nloop_\n"));
		    printf (("_refine_occupancy.class \n"));
		    printf (("_refine_occupancy.details \n"));
		    printf (("_refine_occupancy.treatment \n"));
		    printf (("_refine_occupancy.value \n"));
		    printf (("  ?  ?  ?  ? \n\n"));
		}
	    }
	}
	#
	#   type 3-x remarks

	if ($remark_number ne '  1' && $remark_number ne '  2') {	#???	#???
	    if ($remark_header_flag eq '0') {	#???
		printf (("\nloop_\n"));
		printf (("_database_PDB_remark.id\n"));
		printf (("_database_PDB_remark.text\n"));
		$warning_list{++$warning_flag} = ('#=# DATABASE_PDB_REMARK: '

		  . "Only text in columns 12-70 retained\n");
		++$remark_header_flag;
	    }
	    if ($remark_flag eq '0') {	#???
		printf "%3s\n;", $remark_number;
	    }
	    printf " %-60s \n", $remark_text;
	    ++$remark_flag;
	    ++$flag;
	}
    }

    #===========================================================================
    #  Keyword REVDAT
    #
    #       rev_mod_number	[8-11]	= _database_PDB_rev.num
    #       			= _database_PDB_rev_record.rev_num
    #       rev_cont		[11-12]	= [continuation flag]
    #       rev_date		[14-22]	= _database_PDB_rev.date
    #       			  date in _audit.update_record
    #       rev_name		[24-28]	= _database_PDB_rev_record.details
    #       			  PBD revision name in _audit.update_record
    #       rev_type		[32]	= _database_PDB_rev.mod_type
    #       rev_rec_corr	[40-70]	= _database_PDB_rev_record.type
    #
    #

    if ($first_field eq 'REVDAT') {
	$xlat_flag = $xlat_save;

	$rev_mod_number{$revdat_flag} = substr(($_), 8, 3);
	$rev_date{$revdat_flag} = substr(($_), 14, 9);
	if ($rev_date{$revdat_flag} eq '         ') {
	    $rev_date{$revdat_flag} = $rev_date{$revdat_flag - 1};
	}
	$rev_cont_flag{$revdat_flag} = substr(($_), 11, 2);
	$rev_name{$revdat_flag} = substr(($_), 24, 5);
	$rev_type{$revdat_flag} = substr(($_), 32, 1);
	$rev_rec_corr{$revdat_flag} = substr(($_), 40, 31);
	$rev_mod_of_date{$rev_date{$revdat_flag}} =

	  $rev_mod_number{$revdat_flag} + 0;

	#
	#       The latest revision comes first, take the name
	#       for _audit.revision_id
	#
	if ($revdat_flag == 1) {
	    $aud_rev_id = $rev_name{$revdat_flag};
	    #
	    #       if first REVDAT of a block, save pointer for audit.update_record
	    #
	    ;
	}
	if ($rev_cont_flag{$revdat_flag} eq '  ') {
	    $audit_flag++;
	    $audit_point{$audit_flag} = $revdat_flag;
	    #           prepare OBSLTE and SPRSDE pointer
	    $rev_s_o{$rev_mod_number{$revdat_flag} + 0} = '';
	}

	#   change date format

	{
	    $cal_date = (@dmy = split(/-/, $rev_date{$revdat_flag}, 9999));
	    if ($cal_date >= 1 && $dmy[$cal_date] eq '' && '-' eq ' ') {
		--$cal_date;
	    }
	}
	$rev_date_year{$revdat_flag} = $yyyy{$dmy[3] + 0};
	$rev_date_mon{$revdat_flag} = $mmm2mm{$dmy[2]};
	$rev_date_day{$revdat_flag} = $dmy[1];
	#
	++$revdat_flag;
    }

    #===========================================================================
    #  Keyword SPRSDE and OBSLTE
    #
    #       s_o_type		[1 - 6]	= _database_PDB_rev.status
    #                                     ("obsolete" or ".")
    #       s_o_cont		[9-10]	= [continuation flag]
    #       s_o_date		[12-20]	= _database_PDB_rev.date
    #                               used to find _database_PDB_rev.num
    #       s_o_name		[22-25]	= _database_PDB_rev_record.details
    #       s_o_list		[32-70]	= _database_PDB_rev.replaces (SPRSDE)
    #                                       = _database_PDB_rev.replaced_by (OBSLTE)
    #
    #
    #

    if ($first_field eq 'SPRSDE' || $first_field eq 'OBSLTE') {
	$xlat_flag = $xlat_save;
	$s_o_type{++$s_o_flag} = $first_field;
	$s_o_cont{$s_o_flag} = substr(($_), 9, 2);
	$s_o_date{$s_o_flag} = substr(($_), 12, 9);
	$s_o_name{$s_o_flag} = substr(($_), 22, 4);
	$s_o_list{$s_o_flag} = substr(($_), 32, 39);
    }

    #===========================================================================
    #   Keyword SCALE
    #  CIF provides data items to convert from orthogonal to fractional
    #  coordinates.
    #

    if ($first_field eq 'SCALE1' || $first_field eq 'SCALE2' ||

      $first_field eq 'SCALE3') {
	$xlat_flag = $xlat_save;

	$scale_col1 = substr(($_), 11, 10);
	$scale_col2 = substr(($_), 21, 10);
	$scale_col3 = substr(($_), 31, 10);
	$scale_col4 = substr(($_), 46, 10);
	if ($first_field eq 'SCALE1') {
	    $sc_save{++$scale_flag} = "\n";
	    $sc_save{++$scale_flag} = "loop_\n";
	    $sc_save{++$scale_flag} = "_atom_sites.entry_id\n";
	    $sc_save{++$scale_flag} = "_atom_sites.cartn_transform_axes\n";
	    $sc_save{++$scale_flag} =

	      "_atom_sites.fract_transf_matrix[1][1]\n";
	    $sc_save{++$scale_flag} =

	      "_atom_sites.fract_transf_matrix[1][2]\n";
	    $sc_save{++$scale_flag} =

	      "_atom_sites.fract_transf_matrix[1][3]\n";
	    $sc_save{++$scale_flag} = "_atom_sites.fract_transf_vector[1]\n";
	    $sc_save{++$scale_flag} =

	      "_atom_sites.fract_transf_matrix[2][1]\n";
	    $sc_save{++$scale_flag} =

	      "_atom_sites.fract_transf_matrix[2][2]\n";
	    $sc_save{++$scale_flag} =

	      "_atom_sites.fract_transf_matrix[2][3]\n";
	    $sc_save{++$scale_flag} = "_atom_sites.fract_transf_vector[2]\n";
	    $sc_save{++$scale_flag} =

	      "_atom_sites.fract_transf_matrix[3][1]\n";
	    $sc_save{++$scale_flag} =

	      "_atom_sites.fract_transf_matrix[3][2]\n";
	    $sc_save{++$scale_flag} =

	      "_atom_sites.fract_transf_matrix[3][3]\n";
	    $sc_save{++$scale_flag} = "_atom_sites.fract_transf_vector[3]\n";
	    $sc_save{++$scale_flag} = ('  ' . $head_PDB_code . "\n");
	    $sc_save{++$scale_flag} =

	      "  'See _atom_sites.fract_transf_matrix[i][j]'\n\n";
	}
	$sc_save{++$scale_flag} = ($scale_col1 . ' ' . $scale_col2 . ' ' .

	  $scale_col3 . ' ' . $scale_col4 . "\n");
    }

    #===========================================================================
    #  Keyword SEQADV
    #
    # In the 1995 format, a new record, SEQADV, was added to identify
    # "conflicts between sequence information in the PDB entry and
    # the sequence database" given in DBREF
    #
    #     SEQADV           [1- 6]
    #     seqadv_idcode    [8-11]       = idcode of this entry (not used)
    #     seqadv_resName  [13-15]       = _struct_ref_seq_dif.mon_id
    #     seqadv_chainID   [ 17 ]       = _struct_ref.biol_id
    #     seqadv_seq      [19-23]       = combines seqNum and insertCode
    #                                     used to derive
    #                                     _struct_ref_seq_dif.seq_num 
    #     seqadv_database [25-28]       = _struct_ref.db_name
    #     seqadv_dbAccession
    #                     [30-38]       = _struct_ref.db_code
    #     seqadv_dbRes    [40-42]       = _struct_ref_seq_dif.db_mon_id
    #     seqadv_dbSeq    [44-49]       = _struct_ref_seq_dif.db_seq_num
    #     seqadv_conflict [50-70]       = _struct_ref_seq_dif.details

    if ($first_field eq 'SEQADV') {
	$xlat_flag = $xlat_save;
	$seqadv_resName{++$seqadv_flag} = substr(($_), 13, 3);
	if ($seqadv_resName{$seqadv_flag} eq '   ') {
	    $seqadv_resName{$seqadv_flag} = ' . ';
	}
	$seqadv_chainID{$seqadv_flag} = substr(($_), 17, 1);
	if ($seqadv_chainID{$seqadv_flag} eq ' ') {
	    $seqadv_chainID{$seqadv_flag} = $bcid;
	}
	$seqadv_seq{$seqadv_flag} = substr(($_), 19, 5);
	if ($seqadv_seq{$seqadv_flag} eq '     ') {
	    $seqadv_seq{$seqadv_flag} = '  .  ';
	}
	$seqadv_database{$seqadv_flag} = substr(($_), 25, 4);
	$seqadv_dbAccession{$seqadv_flag} = substr(($_), 30, 9);
	$seqadv_dbRes{$seqadv_flag} = substr(($_), 40, 3);
	$seqadv_dbSeq{$seqadv_flag} = substr(($_), 44, 5);
	$seqadv_conflict{$seqadv_flag} = substr(($_), 50, 21);
	{
	    $numx = (@dblist = split(' ', $seqadv_database{$seqadv_flag},

	      9999));
	    if ($numx >= 1 && $dblist[$numx] eq '' && ' ' eq ' ') {
		--$numx;
	    }
	}
	$seqadv_database{$seqadv_flag} = '.';
	if ($numx > 0) {
	    $seqadv_database{$seqadv_flag} = $dblist[1];
	    for ($j = 2; $j <= $numx; ++$j) {
		$seqadv_database{$seqadv_flag} =

		  ($seqadv_database{$seqadv_flag} . '_' . $dblist[$j]);
	    }
	}
	{
	    $numx = (@dblist = split(' ', $seqadv_dbAccession{$seqadv_flag},

	      9999));
	    if ($numx >= 1 && $dblist[$numx] eq '' && ' ' eq ' ') {
		--$numx;
	    }
	}
	$seqadv_dbAccession{$seqadv_flag} = '.';
	if ($numx > 0) {
	    $seqadv_dbAccession{$seqadv_flag} = $dblist[1];
	    for ($j = 2; $j <= $numx; ++$j) {
		$seqadv_dbAccession{$seqadv_flag} =

		  ($seqadv_dbAccession{$seqadv_flag} . '_' . $dblist[$j]);
	    }
	}
	if ($seqadv_dbRes{$seqadv_flag} eq '   ') {
	    $seqadv_dbRes{$seqadv_flag} = ' . ';
	}
	if ($seqadv_dbSeq{$seqadv_flag} eq '     ') {
	    $seqadv_dbSeq{$seqadv_flag} = '  .  ';
	}
	$seqadv_conflict{$seqadv_flag} = ('Chain ' .

	  $seqadv_chainID{$seqadv_flag} . ': ' .

	  $seqadv_conflict{$seqadv_flag});
	{
	    $numx = (@dblist = split(' ', $seqadv_conflict{$seqadv_flag},

	      9999));
	    if ($numx >= 1 && $dblist[$numx] eq '' && ' ' eq ' ') {
		--$numx;
	    }
	}
	$seqadv_conflict{$seqadv_flag} = '.';
	if ($numx > 0) {
	    $seqadv_conflict{$seqadv_flag} = $dblist[1];
	    for ($j = 2; $j <= $numx; ++$j) {
		$seqadv_conflict{$seqadv_flag} =

		  ($seqadv_conflict{$seqadv_flag} . ' ' . $dblist[$j]);
	    }
	    if ($numx > 1) {
		$seqadv_conflict{$seqadv_flag} = ("'" .

		  $seqadv_conflict{$seqadv_flag} . "'");
	    }
	}
    }

    #==========================================================================
    #  Keyword SEQRES
    #
    #         seq_record_number   [9-10]   = 0 identifies UNK chain
    #         seq_chain_id          [12]   = struct_asym.id 
    #                       used to obtain   _entity_poly_seq.entity_id
    #                       see code in END block
    #         seq_res_count      [14-17]   = count used for UNK chain
    #         seq_text           [20-72]   = _entity_seq_mon_id
    #         seq_number                   = _entity_poly_seq.num
    #         seq_flag                     = num of sequence records
    #

    if ($first_field eq 'SEQRES') {
	$xlat_flag = $xlat_save;

	$seq_record_number{$seqres_flag} = substr(($_), 9, 2);
	$seq_res_count{$seqres_flag} = substr(($_), 14, 4);
	$seq_chain_id{$seqres_flag} = substr(($_), 12, 1);
	$seq_text{$seqres_flag} = substr(($_), 20, 53);

	if ($seq_chain_id{$seqres_flag} eq ' ') {
	    $seq_chain_id{$seqres_flag} = $bcid;
	}
	$seq_entity = $seq_chain_id{$seqres_flag};
	$ent_poly_id{$seq_entity}++;
	if ($ent_poly_id{$seq_entity} == 1) {
	    $next_poly_id = $ent_poly_point{' '};
	    $prev_poly_id = ' ';
	    while ($next_poly_id ne '') {
		$prev_poly_id = $next_poly_id;
		$next_poly_id = $ent_poly_point{$prev_poly_id};
	    }
	    $ent_poly_point{$prev_poly_id} = $seq_entity;
	    $ent_poly_point{$seq_entity} = '';
	    ++$num_poly_ents;
	    $ent_poly_num{$seq_entity} = $num_poly_ents;
	    $entity_seq_num{$seq_entity} = $num_poly_ents;
	    $entities{$num_poly_ents} = $seq_entity;
	}

	++$seqres_flag;
    }

    #=========================================================================
    #   Keyword SHEET
    #
    #
    #   sheet_strand_no        [8-10]  = _struct_sheet_range.id
    #                            _struct_sheet_hbond.range_id_*
    #   sheet_id              [12-14]  = _struct_sheet.id
    #                                    _struct_sheet_hbond.sheet_id
    #                                    _struct_sheet_order.sheet_id
    #                                    _struct_sheet_range.sheet_id
    #   sheet_no_strands      [15-16]  = _struct_sheet.number_strands
    #   sheet_res_name_beg    [18-20]  = _struct_sheet_range.beg_label_comp_id
    #   sheet_chain_id_beg       [22]  = _struct_sheet_range.beg_label_asym_id
    #   sheet_res_seq_beg     [23-27]  = _struct_sheet_range.beg_auth_seq_id
    #   sheet_res_name_end    [29-31]  = _struct_sheet_range.end_label_comp_id
    #   sheet_chain_id_end       [33]  = _struct_sheet_range.end_label_asym_id
    #   sheet_res_seq_end     [34-38]  = _struct_sheet_range.end_auth_seq_id
    #   sheet_sense           [39-40]  = _struct_sheet_order.sense
    #   sheet_atom_name_reg_1 [42-45]  = 
    #                            _struct_sheet_hbond.range_1_beg_label_atom_id
    #   sheet_res_name_reg_1  [46-48]  = 
    #                            _struct_sheet_hbond.range_1_beg_auth_seq_id
    #   sheet_chain_id_reg_1     [50]  = 
    #   sheet_res_seq_reg_1   [51-55]  = 
    #                            _struct_sheet_hbond.range_1_beg_auth_seq_id
    #   sheet_atom_name_reg_2 [57-60]  = 
    #                            _struct_sheet_hbond.range_2_beg_label_atom_id
    #   sheet_res_name_reg_2  [61-63]  = 
    #                            _struct_sheet_hbond.range_2_beg_auth_seq_id
    #   sheet_chain_id_reg_2     [65]  = 
    #   sheet_res_seq_reg_2   [66-70]  = 
    #                            _struct_sheet_hbond.range_2_beg_auth_seq_id
    #
    #   *** note:  The hbond.range_*_end values will be set to the hbond_range*_beg
    #              values, since the PDB format provides only one sample hydrogen
    #              bond for registration

    if ($first_field eq 'SHEET') {
	$xlat_flag = $xlat_save;

	#    Parse field

	$sheet_strand_no{++$sheet_flag} = substr(($_), 8, 3) + 0;
	$sheet_id{$sheet_flag} = substr(($_), 12, 3);
	$sheet_no_strands{$sheet_flag} = substr(($_), 15, 2) + 0;
	$sheet_res_name_beg{$sheet_flag} = substr(($_), 18, 3);
	$sheet_chain_id_beg{$sheet_flag} = substr(($_), 22, 1);
	$sheet_res_seq_beg{$sheet_flag} = substr(($_), 23, 5);
	$sheet_res_name_end{$sheet_flag} = substr(($_), 29, 3);
	$sheet_chain_id_end{$sheet_flag} = substr(($_), 33, 1);
	$sheet_res_seq_end{$sheet_flag} = substr(($_), 34, 5);
	$sheet_sense{$sheet_flag} = substr(($_), 39, 2);
	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = substr(($_), 42, 4);
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $sheet_atom_name_reg_1{$sheet_flag} = $ret_val;
	}
	$sheet_res_name_reg_1{$sheet_flag} = substr(($_), 46, 3);
	$sheet_chain_id_reg_1{$sheet_flag} = substr(($_), 50, 1);
	$sheet_res_seq_reg_1{$sheet_flag} = substr(($_), 51, 5);
	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = substr(($_), 57, 4);
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $sheet_atom_name_reg_2{$sheet_flag} = $ret_val;
	}
	$sheet_res_name_reg_2{$sheet_flag} = substr(($_), 61, 3);
	$sheet_chain_id_reg_2{$sheet_flag} = substr(($_), 65, 1);
	$sheet_res_seq_reg_2{$sheet_flag} = substr(($_), 66, 5);

	if ($sheet_chain_id_beg{$sheet_flag} eq ' ') {
	    $sheet_chain_id_beg{$sheet_flag} = $bcid;
	}
	if ($sheet_chain_id_end{$sheet_flag} eq ' ') {
	    $sheet_chain_id_end{$sheet_flag} = $bcid;
	}
	if ($sheet_sense{$sheet_flag} eq ' 1') {
	    $sheet_sense{$sheet_flag} = 'parallel';
	}
	if ($sheet_sense{$sheet_flag} eq '-1') {
	    $sheet_sense{$sheet_flag} = 'anti-parallel';
	}
	if ($sheet_res_name_reg_1{$sheet_flag} eq '   ' ||

	  $sheet_res_name_reg_1{$sheet_flag} eq '') {
	    $sheet_res_name_reg_1{$sheet_flag} = ' . ';
	}
	if ($sheet_chain_id_reg_1{$sheet_flag} eq ' ') {
	    $sheet_chain_id_reg_1{$sheet_flag} = $bcid;
	}
	if ($sheet_res_seq_reg_1{$sheet_flag} eq '     ' ||

	  $sheet_res_seq_reg_1{$sheet_flag} eq '') {
	    $sheet_res_seq_reg_1{$sheet_flag} = '  .  ';
	}
	if ($sheet_res_name_reg_2{$sheet_flag} eq '   ' ||

	  $sheet_res_name_reg_2{$sheet_flag} eq '') {
	    $sheet_res_name_reg_2{$sheet_flag} = ' . ';
	}
	if ($sheet_chain_id_reg_2{$sheet_flag} eq ' ') {
	    $sheet_chain_id_reg_2{$sheet_flag} = $bcid;
	}
	if ($sheet_res_seq_reg_2{$sheet_flag} eq '     ' ||

	  $sheet_res_seq_reg_2{$sheet_flag} eq '') {
	    $sheet_res_seq_reg_2{$sheet_flag} = '  .  ';
	}
    }

    #=====================================================================
    #   Keyword SIGATM
    #
    #
    #  atom serial number [ 7-11]  = _atom_site.id (used to link to ATOM)
    #  atom type          [13-14]    not used
    #  atom name          [13-16]    not used
    #  alternate location [17]       not used
    #  residue name       [18-20]    not used
    #  chain identifier   [22]       not used
    #  residue seq no.    [23-26]    not used
    #  insertion code     [27]       not used
    #  sd of x-coordinate [31-38]  = _atom_site.cartn_x_esd
    #  sd of y-coordinate [39-46]  = _atom_site.cartn_y_esd
    #  sd of z-coordinate [47-54]  = _atom_site.cartn_z_esd
    #  sd of occupancy    [55-60]  = _atom_site.occupancy_esd
    #  sd of temperature factor 
    #                     [61-66]  = _atom_site.B_iso_or_equiv_esd
    #  footnote number    [68-70]    not used
    #    (February 1992 PDB format)
    #  segment identifier [73-76]    not used
    #    (February 1996 PDB format)
    #  element symbol     [77-78]    not used
    #    (February 1996 PDB format)
    #  charge on atom     [79-80]    not used
    #    (February 1996 PDB format)
    #
    #       
    #

    if ($first_field eq 'SIGATM') {
	$xlat_flag = $xlat_save;

	# parse field

	++$sigatm_flag;

	$sa_atom_serial_number{$sigatm_flag} = substr(($_), 7, 5);
	$sigatm_point{($sa_atom_serial_number{$sigatm_flag} . '|' .

	  $model_flag)} = $sigatm_flag;
	$sig_x{$sigatm_flag} = substr(($_), 31, 8);
	$sig_y{$sigatm_flag} = substr(($_), 39, 8);
	$sig_z{$sigatm_flag} = substr(($_), 47, 8);
	$sig_occ{$sigatm_flag} = substr(($_), 55, 6);
	$sig_temp{$sigatm_flag} = substr(($_), 61, 6);
	if ($sig_x{$sigatm_flag} eq '        ') {
	    $sig_x{$sig_flag} = '    .   ';
	}
	if ($sig_y{$sigatm_flag} eq '        ') {
	    $sig_y{$sigatm_flag} = '    .   ';
	}
	if ($sig_z{$sigatm_flag} eq '        ') {
	    $sig_z{$sigatm_flag} = '    .   ';
	}
	if ($sig_occ{$sigatm_flag} eq '      ') {
	    $sig_occ{$sigatm_flag} = '   .  ';
	}
	if ($sig_temp{$sigatm_flag} eq '      ') {
	    $sig_temp{$sigatm_flag} = '   .  ';
	}
    }

    #=====================================================================
    #   Keyword SIGUIJ
    #
    #
    #       atom serial number   =  matched via pointers to ATOM/HETATM
    #       atom type            =  dropped, taken from ATOM/HETATM
    #       atom name            =  dropped, taken from ATOM/HETATM
    #       alternate location   =  dropped, taken from ATOM/HETATM
    #       residue name         =  dropped, taken from ATOM/HETATM
    #       chain identifier     =  dropped, taken from ATOM/HETATM
    #       residue sequence no. =  dropped, taken from ATOM/HETATM
    #       insertion code       =  dropped, taken from ATOM/HETATM
    #
    #       
    #
    #     Note the different order
    #     PDB           CIF
    # 1.  SIGU[1][1]    SIGU[1][1]
    # 2.  SIGU[2][2]    SIGU[1][2]
    # 3.  SIGU[3][3]    SIGU[1][3]
    # 4.  SIGU[1][2]    SIGU[2][2]
    # 5.  SIGU[1][3]    SIGU[2][3]
    # 6.  SIGU[2][3]    SIGU[3][3]
    #

    if ($first_field eq 'SIGUIJ') {
	$xlat_flag = $xlat_save;

	# parse field

	++$siguij_flag;

	$su_atom_serial_number{$siguij_flag} = substr(($_), 7, 5);
	$siguij_point{($su_atom_serial_number{$siguij_flag} . '|' .

	  $model_flag)} = $siguij_flag;
	$sig_U11{$siguij_flag} = substr(($_), 29, 7);
	$sig_U22{$siguij_flag} = substr(($_), 36, 7);
	$sig_U33{$siguij_flag} = substr(($_), 43, 7);
	$sig_U12{$siguij_flag} = substr(($_), 50, 7);
	$sig_U13{$siguij_flag} = substr(($_), 57, 7);
	$sig_U23{$siguij_flag} = substr(($_), 64, 7);
    }

    #==========================================================================
    #   Keyword SITE
    #
    #       site_seq_no        [8-10]     = _struct_site_gen_id
    #       site_id           [12-14]     = _struct_siite_id
    #                                     = _struct_site_gen_site_id
    #       site_no_res       [16-17]     =
    #       site_res_name_1   [19-21]     = _struct_site_gen_label_comp_id
    #       site_res_name_2   [30-32]
    #       site_res_name_3   [41-43]
    #       site_res_name_4   [52-54]
    #       site_chain_id_1      [23]     = _struct_site_gen_label_asym_id
    #       site_chain_id_2      [34]
    #       site_chain_id_3      [45]
    #       site_chain_id_4      [56]
    #       site_res_seq_1    [24-28]     = _struct_site_gen_auth_seq_id
    #       site_res_seq_2    [35-39]
    #       site_res_seq_3    [46-50]
    #       site_res_seq_4    [57-61]

    if ($first_field eq 'SITE') {
	$xlat_flag = $xlat_save;

	$site_flag_1 = 1;
	$site_seq_no{$site_flag} = substr(($_), 8, 3);
	$site_id{$site_flag} = substr(($_), 12, 3);
	$site_no_res{$site_flag} = substr(($_), 16, 2);
	$site_res_name_1{$site_flag} = substr(($_), 19, 3);
	$site_res_name_2{$site_flag} = substr(($_), 30, 3);
	$site_res_name_3{$site_flag} = substr(($_), 41, 3);
	$site_res_name_4{$site_flag} = substr(($_), 52, 3);
	$site_chain_id_1{$site_flag} = substr(($_), 23, 1);
	$site_chain_id_2{$site_flag} = substr(($_), 34, 1);
	$site_chain_id_3{$site_flag} = substr(($_), 45, 1);
	$site_chain_id_4{$site_flag} = substr(($_), 56, 1);
	$site_res_seq_1{$site_flag} = substr(($_), 24, 5);
	$site_res_seq_2{$site_flag} = substr(($_), 35, 5);
	$site_res_seq_3{$site_flag} = substr(($_), 46, 5);
	$site_res_seq_4{$site_flag} = substr(($_), 57, 5);

	if ($site_chain_id_1{$site_flag} eq ' ') {
	    $site_chain_id_1{$site_flag} = $bcid;
	}
	if ($site_chain_id_2{$site_flag} eq ' ') {
	    $site_chain_id_2{$site_flag} = $bcid;
	}
	if ($site_chain_id_3{$site_flag} eq ' ') {
	    $site_chain_id_3{$site_flag} = $bcid;
	}
	if ($site_chain_id_4{$site_flag} eq ' ') {
	    $site_chain_id_4{$site_flag} = $bcid;
	}
	++$site_flag;
    }

    #===========================================================================
    #   Keyword SLTBRG
    #
    #   Introduced with the February 1996 PDB format
    #
    #        "saltbr"                          = _struct_conn_type.id
    #                                          = _struct_conn.conn_type_id
    #        "defined by user in PDB file"     = _struct_conn_type.criteria
    #        " ? "                             = _struct_conn_type.reference
    #       sb_atom_beg            [13-16]     = _struct_conn.ptnr1_label_atom_id
    #       sb_alt_loc_beg            [17]     = _struct_conn.ptnr1_label_alt_id
    #       sb_res_name_beg        [18-20]     = _struct_conn.ptnr1_label_comp_id
    #       sb_chain_id_beg           [22]     = _struct_conn.ptnr1_label_asym_id
    #       sb_res_seq_num_beg     [23-26]     = _struct_conn.ptnr1_auth_seq_id
    #       sb_icode_beg              [27]     append to
    #                                            _struct_conn.ptnr1_auth_seq_id
    #       sb_atom_end            [43-46]     = _struct_conn.ptnr2_label_atom_id
    #       sb_alt_loc_end            [47]     = _struct_conn.ptnr2_label_alt_id
    #       sb_res_name_end        [48-50]     = _struct_conn.ptnr2_label_comp_id
    #       sb_chain_id_end           [52]     = _struct_conn.ptnr2_label_asym_id
    #       sb_res_seq_num_end     [53-56]     = _struct_conn.ptnr2_auth_seq_id
    #       sb_icode_end              [57]     append to
    #                                            _struct_conn.ptnr2_auth_seq_id
    #       sb_symop1              [60-65]     = _struct_conn.ptnr1_symmetry
    #       sb_symop2              [67-72]     = _struct_conn.ptnr2_symmetry
    #

    if ($first_field eq 'SLTBRG') {
	$xlat_flag = $xlat_save;

	$sb_atom_beg{++$sltbrg_flag} = substr(($_), 13, 4);
	$sb_alt_loc_beg{$sltbrg_flag} = substr(($_), 17, 1);
	$sb_res_name_beg{$sltbrg_flag} = substr(($_), 18, 3);
	$sb_chain_id_beg{$sltbrg_flag} = substr(($_), 22, 1);
	$sb_res_seq_num_beg{$sltbrg_flag} = substr(($_), 23, 5);
	$sb_atom_end{$sltbrg_flag} = substr(($_), 43, 4);
	$sb_alt_loc_end{$sltbrg_flag} = substr(($_), 47, 1);
	$sb_res_name_end{$sltbrg_flag} = substr(($_), 48, 3);
	$sb_chain_id_end{$sltbrg_flag} = substr(($_), 52, 1);
	$sb_res_seq_num_end{$sltbrg_flag} = substr(($_), 53, 5);
	$sb_symm_1{$sltbrg_flag} = substr(($_), 60, 6);
	$sb_symm_2{$sltbrg_flag} = substr(($_), 67, 6);

	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = $sb_atom_beg{$sltbrg_flag};
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $temp_name = $ret_val;
	}
	$sb_atom_beg{$sltbrg_flag} = $temp_name;
	{
	    #
	    #
	    # fix up atom_name by squeezing out blanks in the middle
	    #
	    $temp_a_name = $sb_atom_end{$sltbrg_flag};
	    if (substr($temp_a_name, 3, 1) eq ' ') {
		$temp_a_name = (substr($temp_a_name, 1,

		  2) . substr($temp_a_name, 4, 1) . ' ');
	    }
	    if (substr($temp_a_name, 2, 1) eq ' ') {
		$temp_a_name = (' ' . substr($temp_a_name, 1,

		  1) . substr($temp_a_name, 3, 2));
	    }
	    if ($temp_a_name eq '    ') {
		$temp_a_name = ' .  ';
	    }
	    $ret_val = $temp_a_name;

	    $temp_name = $ret_val;
	}
	$sb_atom_end{$sltbrg_flag} = $temp_name;
	if ($sb_alt_loc_beg{$sltbrg_flag} eq ' ') {
	    $sb_alt_loc_beg{$sltbrg_flag} = '.';
	}
	if ($sb_alt_loc_end{$sltbrg_flag} eq ' ') {
	    $sb_alt_loc_end{$sltbrg_flag} = '.';
	}
	if ($sb_res_name_beg{$sltbrg_flag} eq '   ') {
	    $sb_res_name_beg{$sltbrg_flag} = ' . ';
	}
	if ($sb_res_name_end{$sltbrg_flag} eq '   ') {
	    $sb_res_name_end{$sltbrg_flag} = ' . ';
	}
	if ($sb_chain_id_beg{$sltbrg_flag} eq ' ') {
	    $sb_chain_id_beg{$sltbrg_flag} = $bcid;
	}
	if ($sb_chain_id_end{$sltbrg_flag} eq ' ') {
	    $sb_chain_id_end{$sltbrg_flag} = $bcid;
	}
	if ($sb_symm_1{$sltbrg_flag} eq '      ') {
	    $sb_symm_1{$sltbrg_flag} = '   .   ';
	}
	else {
	    $sb_symm_1{$sltbrg_flag} = (substr($sb_symm_1{$sltbrg_flag}, 1,

	      3) . '_' . substr($sb_symm_1{$sltbrg_flag}, 4, 3));
	}
	if ($sb_symm_2{$sltbrg_flag} eq '      ') {
	    $sb_symm_2{$sltbrg_flag} = '   .   ';
	}
	else {
	    $sb_symm_2{$sltbrg_flag} = (substr($sb_symm_2{$sltbrg_flag}, 1,

	      3) . '_' . substr($sb_symm_2{$sltbrg_flag}, 4, 3));
	}
    }

    #===========================================================================
    #  Keyword SOURCE
    #
    # The PDB describes the source of all components of the structure here
    # and has no way to parse out the individual entities in the 1992 format
    # the situation changes with a use of keywords in the 1995 format
    # all of this is added to _struct.title (see code after COMPND)
    #
    #     source_cont      [9-10]
    #     source_text     [11-70]       = _struct.title

    if ($first_field eq 'SOURCE') {
	$xlat_flag = $xlat_save;

	$source_cont = substr(($_), 9, 2);
	$source_text = substr(($_), 11, 60);
	if ($convtext eq 'yes') {
	    #
	    # apply PDB typsetting codes if any to a line
	    #
	    {
		$lx_tl = length($source_text);
		$tx_tl = $source_text;
		$lostr = '';
		for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
		    $cx_tl = substr($tx_tl, $ix_tl, 1);
		    $cx_tl = substr(($lcaz . $cx_tl), index(($UCAZ . $cx_tl),

		      $cx_tl), 1);
		    $lostr = ($lostr . $cx_tl);
		}
	    }

	    $lstr = length($lostr);
	    $mystr = '';
	    $pchar = ' ';
	    for ($qtsi = 1; $qtsi <= $lstr; ++$qtsi) {
		$mychar = substr($lostr, $qtsi, 1);
		if ($pchar eq ' ' || $pchar eq ',' || $pchar eq '.' ||

		  $pchar eq '(' || $pchar eq '*' || $pchar eq '/') {
		    $mychar = substr(($UCAZ . $mychar),

		      index(($lcaz . $mychar), $mychar), 1);
		}
		if (($mychar ne '*' && $mychar ne "\$" && $mychar ne '/') ||

		  ($mychar eq $pchar)) {	#???
		    $mystr = ($mystr . $mychar);
		}
		if ($pchar eq '/') {
		    if ($mychar eq "\$" || $mychar eq '-') {
			$pchar = $mychar;
		    }
		}
		else {
		    $pchar = $mychar;
		}
	    }
	    $ret_val = $mystr;

	    $source_text = $ret_val;
	}

	if ($source_cont eq '  ') {
	    $compnd{$compnd_flag++} = 'Source::';
	}
	$bp = '     ';
	if ($source_cont ne '  ') {
	    $bp = '    ';
	}
	$compnd{$compnd_flag++} = ($bp . $source_text);
    }

    #===========================================================================
    #   Keyword SSBOND
    #
    #   Note, in the 1996 PDB format description, the former comment field
    #   in columns 41-70 was removed and columns 60-65 and 67-72 used for
    #   symmetry operators.  It is believed that the comment field was not
    #   used for any significant number of entries in the past.  Therefore
    #   the symmetry operators will be parsed for all entries.
    #
    #        "disulf"                          = _struct_conn_type.id
    #                                          = _struct_conn.conn_type_id
    #        "defined by user in PDB file"     = _struct_conn_type.criteria
    #        " ? "                             = _struct_conn_type.reference
    #       ssflag                             = _struct_conn.id
    #       ssbond_num              [8-10]     =
    #       ssbond_res_name_beg    [12-14]     = _struct_conn.ptnr1_label_comp_id
    #       ssbond_chain_id_beg       [16]     = _struct_conn.ptnr1_label_asym_id
    #       ssbond_res_seq_num_beg [18-22]     = _struct_conn.ptnr1_auth_seq_id
    #       	n/a			   = _struct_conn.ptnr1_label_alt_id
    #       ssbond_res_name_end    [26-28]     = _struct_conn.ptnr2_label_comp_id
    #       ssbond_chain_id_end       [30]     = _struct_conn.ptnr2_label_asym_id
    #       ssbond_res_seq_num_end [32-36]     = _struct_conn.ptnr2_auth_seq_id
    #               n/a                        = _struct_conn.ptnr2_label_alt_id
    #       ssbond_symm_1          [60-65]     = _struct_conn.ptnr1_symmetry
    #       ssbond_symm_2          [67-72]     = _struct_conn.ptnr2_symmetry
    #

    if ($first_field eq 'SSBOND') {
	$xlat_flag = $xlat_save;

	$ssbond_num{$ssbond_flag} = substr(($_), 8, 3);
	$ssbond_res_name_beg{$ssbond_flag} = substr(($_), 12, 3);
	$ssbond_chain_id_beg{$ssbond_flag} = substr(($_), 16, 1);
	$ssbond_res_seq_num_beg{$ssbond_flag} = substr(($_), 18, 5);
	$ssbond_res_name_end{$ssbond_flag} = substr(($_), 26, 3);
	$ssbond_chain_id_end{$ssbond_flag} = substr(($_), 30, 1);
	$ssbond_res_seq_num_end{$ssbond_flag} = substr(($_), 32, 5);
	$ssbond_symm_1{$ssbond_flag} = substr(($_), 60, 6);
	$ssbond_symm_2{$ssbond_flag} = substr(($_), 67, 6);

	if ($ssbond_chain_id_beg{$ssbond_flag} eq ' ') {
	    $ssbond_chain_id_beg{$ssbond_flag} = $bcid;
	}
	if ($ssbond_chain_id_end{$ssbond_flag} eq ' ') {
	    $ssbond_chain_id_end{$ssbond_flag} = $bcid;
	}
	if ($ssbond_symm_1{$ssbond_flag} eq '      ') {
	    $ssbond_symm_1{$ssbond_flag} = '   .   ';
	}
	else {
	    $ssbond_symm_1{$ssbond_flag} =

	      (substr($ssbond_symm_1{$ssbond_flag}, 1,

	      3) . '_' . substr($ssbond_symm_1{$ssbond_flag}, 4, 3));
	}
	if ($ssbond_symm_2{$ssbond_flag} eq '      ') {
	    $ssbond_symm_2{$ssbond_flag} = '   .   ';
	}
	else {
	    $ssbond_symm_2{$ssbond_flag} =

	      (substr($ssbond_symm_2{$ssbond_flag}, 1,

	      3) . '_' . substr($ssbond_symm_2{$ssbond_flag}, 4, 3));
	}

	++$ssbond_flag;
    }

    #=======================================================================
    #   Keyword TER
    #
    #  Used here to increment the number TER records
    #  See ATOM/HETATM for rest of processing.
    #
    #       ter_num         [7-11]     =
    #       ter_res_name   [18-20]     =
    #       ter_chain_id      [22]     =
    #       ter_res_seq_num[23-26]     =

    if ($first_field eq 'TER') {
	$ter_flag++;
    }

    #========================================================================
    #  Keyword TITLE
    #
    # In the 1995 format, the record type TITLE was added, in many ways
    # replacing the former fucntion of COMPND, which became more detailed
    # with a keyword structure.
    #
    #     title_cont      [9-10]
    #     title_text     [11-70]       = _struct.title

    if ($first_field eq 'TITLE') {
	$xlat_flag = $xlat_save;

	$title_cont = substr(($_), 9, 2);
	$title_text = substr(($_), 11, 60);
	$bp = '     ';
	if ($title_cont ne '  ') {
	    $bp = '    ';
	}
	$compnd{$compnd_flag++} = ($bp . $title_text);
    }

    #=======================================================================
    #   Keyword TURN
    #
    # From the PDB file it is not possible to determine _struct_conf.conf_type_id
    #
    # * indicates a _struct_topol data item which is not currently used
    #
    #            "TURN"                = _struct_conf.conf_type_id &
    #           unknown                = _struct_conf_conf_type_id
    #          "From PDB"              = _struct_conf_type_criteria
    #             "?"                  = _struct_conf_type_reference
    #           "Turn"                 = _struct_topol_type *
    #          "From PDB"              = _struct_topol_criteria *
    #      turn_num             [8-10] = _struct_conf.id
    #      turn_id             [12-14] = _struct_topol_id *
    #      turn_res_name_beg   [16-18] = _struct_conf.beg_label_comp_id
    #      turn_chain_id_beg      [20] = _struct_conf.beg_label_asym_id
    #      turn_res_seq_beg [21-24] = _struct_conf.beg_auth_seq_id
    #      turn_res_name_end   [27-29] = _struct_conf.end_label_comp_id
    #      turn_chain_id_end      [31] = _struct_conf.end_label_asym_id
    #      turn_res_seq_end [32-35] = _struct_conf.end_auth_seq_id
    #      turn_comment        [41-70] = _struct_conf.details

    if ($first_field eq 'TURN') {
	$xlat_flag = $xlat_save;

	# parse field

	$turn_num{$ss_flag} = substr(($_), 8, 3);
	$turn_id{$ss_flag} = substr(($_), 12, 3);
	$turn_res_name_beg{$ss_flag} = substr(($_), 16, 3);
	$turn_chain_id_beg{$ss_flag} = substr(($_), 20, 1);

	$turn_res_seq_beg{$ss_flag} = substr(($_), 21, 5);
	$turn_res_name_end{$ss_flag} = substr(($_), 27, 3);
	$turn_chain_id_end{$ss_flag} = substr(($_), 31, 1);

	$turn_res_seq_end{$ss_flag} = substr(($_), 32, 5);
	$turn_comment{$ss_flag} = substr(($_), 41, 30);

	# strip blanks from id
	{
	    $num_x = (@xxx = split(' ', $turn_id{$ss_flag}, 9999));
	    if ($num_x >= 1 && $xxx[$num_x] eq '' && ' ' eq ' ') {
		--$num_x;
	    }
	}
	$turn_id{$ss_flag} = '';
	if ($num_x == 1) {
	    $turn_id{$ss_flag} = $xxx[1];
	}
	if ($num_x == 2) {
	    $turn_id{$ss_flag} = ($xxx[1] . '_' . $xxx[2]);
	}
	if ($turn_chain_id_beg{$ss_flag} eq ' ') {
	    $turn_chain_id_beg{$ss_flag} = $bcid;
	}
	if ($turn_chain_id_end{$ss_flag} eq ' ') {
	    $turn_chain_id_end{$ss_flag} = $bcid;

	    # give real name to turn classes (should all be TURN_P, but might as
	    # well check)

	    ;
	}
	$t_class_suffix = '_P';
	{
	    $num_x = (@xxx = split(' ',

	      ($turn_res_name_end{$ss_flag} . ' ' .

	      $turn_res_name_beg{$ss_flag}), 9999));
	    if ($num_x >= 1 && $xxx[$num_x] eq '' && ' ' eq ' ') {
		--$num_x;
	    }
	}
	foreach $i ($[ .. $#na_list) {
	    if ($na_list[$i] eq $xxx[1] || $na_list[$i] eq $xxx[2]) {	#???	#???
		$t_class_suffix = '_N';
	    }
	}

	$turn_class{$ss_flag} = ('TURN' . $t_class_suffix);

	++$ss_flag;
	++$turn_flag;
    }

    #============================================================================
    #  Keyword TVECT
    #
    # _database_pdb_tvect.id
    # _database_pdb_tvect.vector[1]  .. _3
    # _databsee_pdb_tvect.details

    if ($first_field eq 'TVECT') {
	$xlat_flag = $xlat_save;

	$tvect_id = substr(($_), 8, 3);
	$tvect_ent1 = substr(($_), 11, 10);
	$tvect_ent2 = substr(($_), 21, 10);
	$tvect_ent3 = substr(($_), 31, 10);
	$tvect_com = substr(($_), 41, 30);

	# print loop headers

	if ($tvect_flag eq '0') {	#???
	    $tv_save{++$tvect_flag} = "\n\n\n";

	    $tv_save{++$tvect_flag} = "\nloop_ \n";
	    $tv_save{++$tvect_flag} = "_database_pdb_tvect.id\n";
	    $tv_save{++$tvect_flag} = "_database_pdb_tvect.vector[1] \n";
	    $tv_save{++$tvect_flag} = "_database_pdb_tvect.vector[2] \n";
	    $tv_save{++$tvect_flag} = "_database_pdb_tvect.vector[3] \n";
	    $tv_save{++$tvect_flag} = "_database_pdb_tvect.details \n";
	}
	if ($tvect_id eq '   ') {
	    $tvect_id = '  0';
	}
	$tv_save{++$tvect_flag} = (' ' . $tvect_id . ' ' . $tvect_ent1 . ' ' .

	  $tvect_ent2 . ' ' . $tvect_ent3 . ' ' . "\n; " . $tvect_com .

	  "\n;\n");
    }
}

#  check seqres info for duplicates
for ($i = 1; $i <= $num_poly_ents; ++$i) {
    $seq_chain{$i} = '';
    $seq_sig{$i} = '0000000000';
}
for ($i = 1; $i <= $num_poly_ents; ++$i) {
    $entities_list{$i} = ('Chain: ' . $entities{$i});
}
for ($is = 1; $is < $seqres_flag; ++$is) {
    {
	$num_seq = (@seq_res = split(' ', $seq_text{$is}, 9999));
	if ($num_seq >= 1 && $seq_res[$num_seq] eq '' && ' ' eq ' ') {
	    --$num_seq;
	}
    }
    $stutter = 0;
    if ($seq_record_number{$is} eq ' 0') {
	$num_seq = $seq_res_count{$is};
	$stutter = 1;
    }
    for ($i = 1; $i <= $num_seq; ++$i) {
	$k = $ent_poly_num{$seq_chain_id{$is}};
	$isx = $i * (1 - $stutter) + $stutter;
	$cur_res = $seq_res[$isx];
	$x_cur_res = substr(('   ' . $cur_res), 1 + length($cur_res), 3);
	++$res_res_flag;
	if ($res_count{$x_cur_res} eq '') {
	    $res_count{$x_cur_res} = 0;
	    $res_formul{$x_cur_res} = '.';
	    $res_name{$x_cur_res} = '.';
	    ++$num_OTRESES;
	    $res_id{$num_AARESES + $num_NARESES + $num_OTRESES} = $x_cur_res;
	}
	++$res_count{$x_cur_res};
	$found = 'no';
	foreach $j ($[ .. $#na_list) {
	    if ($na_list[$j] eq $cur_res) {	#???
		if (substr($entities_list{$k}, 1, 3) eq 'Pro') {
		    $entities_list{$k} = ('Protein/Nucleic Acid chain: ' .

		      $entities{$k});
		}
		else {
		    $entities_list{$k} = ('Nucleic Acid chain: ' .

		      $entities{$k});
		}
		$found = 'yes';
		last;
	    }
	}
	if ($found eq 'no') {
	    foreach $j ($[ .. $#aa_list) {
		if ($aa_list[$j] eq $cur_res) {	#???
		    if (index($entities_list{$k}, 'Nuc') > 0) {
			$entities_list{$k} = ('Protein/Nucleic Acid chain: ' .

			  $entities{$k});
		    }
		    else {
			$entities_list{$k} = ('Protein chain: ' .

			  $entities{$k});
		    }
		    last;
		}
	    }
	}
	if ($res_code{$cur_res} eq '') {
	    ++$num_res_name;
	    if ($num_res_name > length($charx)) {
		$num_res_name -= length($charx);
	    }
	    $res_code{$cur_res} = substr($charx, $num_res_name, 1);
	}
	$seq_chain{$ent_poly_num{$seq_chain_id{$is}}} =

	  ($seq_chain{$ent_poly_num{$seq_chain_id{$is}}} .

	  $res_code{$cur_res});
    }
}
for ($i = 1; $i <= $num_poly_ents; ++$i) {
    for ($j = 1; $j < length($seq_chain{$i}); ++$j) {
	$code_pair = substr($seq_chain{$i}, $j, 2);
	if ($pair_code{$code_pair} eq '') {
	    ++$num_res_pair;
	    if ($num_res_pair > length($numl)) {
		$num_res_pair -= length($numl);
	    }
	    $pair_code{$code_pair} = substr($numl, $num_res_pair, 1);
	}
	$seq_dig = substr($seq_sig{$i}, $pair_code{$code_pair} + 1, 1);
	$seq_dig = $seq_dig + 1;
	if ($seq_dig > 9) {
	    $seq_dig -= 10;
	}
	$seq_sig{$i} = (substr($seq_sig{$i}, 1,

	  $pair_code{$code_pair}) . $seq_dig . substr($seq_sig{$i},

	  $pair_code{$code_pair} + 2,

	  length($seq_sig{$i}) - $pair_code{$code_pair} - 1));
    }
}
for ($i = 1; $i <= $num_poly_ents; ++$i) {
    $el_point{$i} = length($entities_list{$i});
}
$known_links = '';
for ($i = 1; $i < $num_poly_ents; ++$i) {
    for ($j = $i + 1; $j <= $num_poly_ents; ++$j) {
	{
	    $lsal = length($seq_sig{$i});
	    $seqlim = $lsal / 2;
	    if ($seqlim > 6) {
		$seqlim = 6;
	    }
	    $lssr = length($seq_sig{$j});
	    $isr = 1;
	    $nm = 0;
	    $missed = 1;
	    for ($isha = 1; $isha <= $lsal; ++$isha) {
		$seq_match{$isha} = 0;
	    }
	    for ($ipass = 0; (($ipass < 4) && ($missed != 0)); $ipass++) {
		$missed = 0;
		$isha = 1;
		$isr = 1;
		while ($isha <= $lsal) {	#???
		    $sm = 0;
		    $lm = $lsal - $isha + 1;
		    if ($seq_match{$isha} == 0) {
			while ($lm > 0 && $sm == 0) {
			    $cr = substr($seq_sig{$i}, $isha, $lm);
			    $sm = index(substr($seq_sig{$j}, $isr,

			      $lssr - $isr + 1), $cr);
			    $sf = 1;
			    if (!($known_links eq '') && !($sm == 0) &&

			      ($lm < $seqlim)) {	#???
				if ($isha < $lsal && $ipass > 0) {	#???
				    $kl = substr($known_links,

				      $isha + $lm - 1, 1);
				    if ($kl eq '1' && $ipass < 3 &&

				      $seq_match{$isha + $lm} > 0) {
					if (!($sm + $isr + $lm - 1 ==

					  $seq_match{$isha + $lm}) &&

					  ($ipass < 3)) {
					    $sm = 0;
					}
				    }
				    else {
					$kbase = $isha + $lm - 1;
					$kbreak = 0;
					while (($kbase < $lsal) &&	#???

					  !($seq_match{$kbase + 1} > 0)) {
					    $kbase++;
					    if (substr($known_links, $kbase,

					      1) eq '0') {
						$kbreak++;
					    }
					}
					if (($seq_match{$kbase + 1} > 0) &&

					  ((($sm + $isr - 1 + $kbase - $isha -

					  $lm > $seq_match{$kbase + 1} +

					  $kbreak * 3 - 1)) ||

					  ($sm + $isr - 1 + $kbase - $isha -

					  $lm < $seq_match{$kbase + 1} -

					  $kbreak * 3 - 2 - $seqlim * 2))) {
					    $sf = -1;
					}
				    }
				}
				if ($isha > 1 && !($sm == 0) && $sf == 1) {
				    $kl = substr($known_links, $isha, 1);
				    if ($kl eq '1' && $ipass < 3 &&

				      $seq_match{$isha - 1} > 0) {
					if (!($sm + $isr - 2 ==

					  $seq_match{$isha - 1}) &&

					  ($ipass < 3)) {
					    $sm = 0;
					}
				    }
				    else {
					$kbase = $isha;
					$kbreak = 0;
					while ((($kbase - 1) > 0) &&

					  !($seq_match{$kbase} > 0)) {
					    $kbase--;
					    if (substr($known_links, $kbase,

					      1) eq '0') {
						$kbreak++;
					    }
					}
					if (($seq_match{$kbase} > 0) &&

					  ((($sm + $isr - 1 + $kbase - $isha <

					  $seq_match{$kbase} - $kbreak * 3 -

					  1)) || $sm + $isr - 1 + $kbase -

					  $isha > $seq_match{$kbase} + $kbreak

					  * 3 + $seqlim * 2 - 1)) {
					    $sf = -1;
					}
				    }
				}
			    }

			    if ($sm > 0) {
				$kbase = $isha;
				while (($kbase > 1) &&

				  ($seq_match{$kbase - 1} <= 0) &&

				  ($sf > 0) &&

				  (substr($seq_sig{$i}, $kbase - 1,

				  1) eq substr($seq_sig{$j},	#???

				  $sm + $isr - 2 + $kbase - $isha, 1))) {
				    $kbase--;
				    $nm++;
				    $seq_match{$kbase} = $sm + $isr - 1 +

				      $kbase - $isha;
				}
				for ($ksa = 0; $ksa < length($cr); ++$ksa) {
				    if ($seq_match{$isha} == 0) {
					$seq_match{$isha} = ($sm + $ksa + $isr

					  - 1) * $sf;
					++$isha;
				    }
				    else {
					$cr = substr($cr, 1, $ksa);
				    }
				}
				$isr = $sm + length($cr) + $isr - 1;
				if (length($cr) > 3 || $isha == $lsal + 1) {
				    if ($sf > 0) {
					$nm += length($cr);
				    }
				}
				last;
			    }
			    else {
				if ($lm > 16) {
				    $lm = int($lm * .707);
				}
				else {
				    --$lm;
				    if ($ipass < 3 && $lm < $seqlim) {	#???
					if (($isha == 1) ||

					  ($isha > 1 &&

					  $seq_match{$isha - 1} == 0) &&

					  ($isha < $lsal - $lm + 1 &&

					  $seq_match{$isha + $lm} == 0)) {
					    $lm = 0;
					}
				    }
				}
			    }
			}
		    }
		    if ($sm == 0) {
			++$isha;
			$missed++;
		    }
		}
	    }
	    $qxmaxlen = $lsal;
	    if ($lssr > $lsal) {	#???
		$qxmaxlen = $lssr;
	    }
	    $ret_val = 100 * ($nm / $qxmaxlen);

	    $xcomp = $ret_val;
	}

	if ($seq_sig{$i} eq $seq_sig{$j}) {	#???
	    if ($seq_chain{$i} eq $seq_chain{$j}) {	#???
		$entity_seq_num{$entities{$j}} =

		  $entity_seq_num{$entities{$i}};
		$el_point{$i} += 3;
		$xtemp = ', ';
		if ($el_point{$i} >= 70) {
		    $xtemp = ",\n  ";
		    $el_point{$i} = 3;
		}
		$entities_list{$i} = ($entities_list{$i} . $xtemp .

		  $entities{$j});
	    }
	}
	else {
	    if ($xcomp > 85) {
		printf

		  "# **** WARNING **** approx %5d%% homology %3s to %3s \n",

		  int($xcomp), $entities{$i}, $entities{$j};
		$warning_list{++$warning_flag} =

		  sprintf("#=# ENTITY_POLY_SEQ: Approx %5d%%%% homology %3s to %3s \n",

		  int($xcomp), $entities{$i}, $entities{$j});
	    }
	}
    }
}
#   build a complete _entity list from ATOM HETATM and FORMUL records
$num_ents = $num_poly_ents + 1;
$next_non_poly = ' ';
$next_non_poly = $ent_non_poly_point{$next_non_poly};
while ($next_non_poly ne '') {
    if ($entity_seq_num{$next_non_poly} eq '' ||

      $entity_seq_num{$next_non_poly} + 0 > $num_poly_ents) {
	$entities{$num_ents} = $next_non_poly;
	$entity_seq_num{$next_non_poly} = $num_ents;
	++$num_ents;
    }
    $next_non_poly = $ent_non_poly_point{$next_non_poly};
}
{
    # Process AUDIT information
    $| = 1;
    printf (("\n\n\n"));
    printf (("####################\n"));
    printf (("#                  #\n"));
    printf (("# AUDIT            #\n"));
    printf (("#                  #\n"));
    printf (("####################\n\n\n"));
    printf "_audit.revision_id         %4s\n", $aud_rev_id;
    {
	$null = (@dmy = split(/-/, $head_dep_date, 9999));
	if ($null >= 1 && $dmy[$null] eq '' && '-' eq ' ') {
	    --$null;
	}
    }
    printf "_audit.creation_date       %4s-%2s-%2s\n", $yyyy{$dmy[3] + 0},

      $mmm2mm{$dmy[2]}, $dmy[1];
    printf (("_audit.update_record\n; "));
    if ($audit_flag > 1) {
	for ($j = 1; $j < $audit_flag; ++$j) {
	    $i = $audit_point{$audit_flag - $j};
	    printf "%4s-%2s-%2s  PDB revision %5s\n  ", $rev_date_year{$i},

	      $rev_date_mon{$i}, $rev_date_day{$i}, $rev_name{$i};
	}
	# end for (j=1; j < audit_flag; ++j)

	;
    }
    # end if (audit_flag > 1 )
    $x_time = localtime();
    if ($x_time eq '') {	#???
	$x_time =

	  system(("/bin/echo `date +%Y-%m-%d` ' Converted to mmCIF format by pdb2cif.pl'"

	  . ' ' . '2.4.2'));
    }
    else {
	printf (((substr($x_time, 21, 4) . '-' . $mmm2mm{substr($x_time, 5,

	  3)} . '-' . substr(substr($x_time, 9, 2) + 100, 2,

	  2) . ('  Converted to mmCIF format by pdb2cif.pl' . ' ' .

	  "2.4.2\n"))));
    }
    # If this call to system does not work on your system
    # First try changing the "+%Y" to "+19%y" and if that fails, just         #
    # comment out the two lines and forget it.  All they do is add a trace    #
    # of when the conversion  was done                                        #
    ###########################################################################
    $| = 0;
    if ($xlat_flag > 0) {
	printf (("***** Did Not Translate:\n"));
	for ($X = 1; $X <= $xlat_flag; ++$X) {
	    printf ((($non_xlated{$X} . "\n")));
	}
	$warning_list{++$warning_flag} =

	  sprintf("#=# AUDIT: Did not translate %6d PDB records\n",

	  $xlat_flag);
    }

    printf ((';'));
}
for ($i = 1; $i <= $keywrd_flag; ++$i) {
    printf (($key_save{$i}));
}
# Print out _struct_biol_info

if ($verbose eq 'yes') {
    printf (("\nloop_\n"));
    printf (("_struct_biol_keywords.biol_id\n"));
    printf (("_struct_biol_keywords.text\n"));
    printf " %4s   ?\n", $head_PDB_code;

    printf (("\nloop_\n"));
    printf (("_struct_biol_view.biol_id\n"));
    printf (("_struct_biol_view.id\n"));
    printf (("_struct_biol_view.details\n"));
    printf (("_struct_biol_view.rot_matrix[1][1]\n"));
    printf (("_struct_biol_view.rot_matrix[1][2]\n"));
    printf (("_struct_biol_view.rot_matrix[1][3]\n"));
    printf (("_struct_biol_view.rot_matrix[2][1]\n"));
    printf (("_struct_biol_view.rot_matrix[2][2]\n"));
    printf (("_struct_biol_view.rot_matrix[2][3]\n"));
    printf (("_struct_biol_view.rot_matrix[3][1]\n"));
    printf (("_struct_biol_view.rot_matrix[3][2]\n"));
    printf (("_struct_biol_view.rot_matrix[3][3]\n\n"));
    printf " %4s 1  ? ? ? ? ? ? ? ? ? ?\n", $head_PDB_code;
}
# end if (verbose == "yes")

# Information about entities which are polymers
#
#  seqres_flag       count of seqres records + 1
#  is                index from 1 to seqres_flag-1
#  seq_chain_id[is]  chainID from currect SEQRES record (e.g. A, B, *)
#  entity_seq_num[chainID]
#                    the entity number for chainID
#  ent_poly_num[chainID]
#                    the PDB component number of the sequence
#                    this number counts chains, not entities
#  seq_start_num[entityNum]
#                    the value of _entity_poly_seq.num at which the
#                    sequence for entityNun begins
#  seq_end_num[entityNum]
#                    the value of _entity_poly_seq.num at which the
#                    sequence for entityNun ends
#  seq_number        counter for _entity_poly_seq.num
#  entity_poly_seq[seqNum]
#                    value of _entity_poly_seq.mon_id for
#                    _entity_poly_seq.num = seqNum

if ($seqres_flag + 0 > 1) {
    printf (("\n"));
    printf (("##########################\n"));
    printf (("#                        #\n"));
    printf (("# ENTITY_POLY_SEQ        #\n"));
    printf (("#                        #\n"));
    printf (("##########################\n\n"));
    printf (("loop_ \n"));
    printf (("_entity_poly_seq.entity_id\n"));
    printf (("_entity_poly_seq.num\n"));
    printf (("_entity_poly_seq.mon_id\n"));
    $iscount = 0;
    $iscurr = $seq_chain_id{1};
    $seq_number = 1;
    $seq_start_num{$entity_seq_num{$iscurr}} = 1;
    $seq_end_num{$entity_seq_num{$iscurr}} = 1;
    for ($is = 1; $is < $seqres_flag; ++$is) {
	if ($ent_poly_num{$seq_chain_id{$is}} eq	#???

	  $entity_seq_num{$seq_chain_id{$is}}) {
	    if ($iscurr ne $seq_chain_id{$is}) {	#???
		$seq_start_num{$entity_seq_num{$seq_chain_id{$is}}} =

		  $seq_number;
		if ($iscount > 0) {
		    printf (("\n"));
		    $iscount = 0;
		}
	    }
	    $iscurr = $seq_chain_id{$is};
	    {
		$num_seq = (@seq_res = split(' ', $seq_text{$is}, 9999));
		if ($num_seq >= 1 && $seq_res[$num_seq] eq '' && ' ' eq ' ') {
		    --$num_seq;
		}
	    }
	    $stutter = 0;
	    if ($seq_record_number{$is} eq ' 0') {
		$num_seq = $seq_res_count{$is};
		$stutter = 1;
	    }
	    for ($i = 1; $i <= $num_seq; ++$i) {
		$isx = $i * (1 - $stutter) + $stutter;
		printf ' %5d %4s %3s', $ent_poly_num{$seq_chain_id{$is}},

		  $seq_number, $seq_res[$isx];
		$entity_poly_seq{$seq_number} = $seq_res[$isx];
		$seq_end_num{$entity_seq_num{$iscurr}} = $seq_number;
		++$seq_number;
		++$iscount;
		if ($iscount > 4) {
		    printf (("\n"));
		    $iscount = 0;
		}
	    }
	}
    }
    if ($iscount > 0) {
	printf (("\n"));
    }
}

#
#   process _entity info combining DBREF SEQADV ATOM and HETATM records
#
$oxt_found = 0;
$poxt_found = 0;
$prior_res_name = '';
$pprior_res_name = '';
$prior_block = '';
$pprior_block = '';
$prior_x = 0;
$prior_y = 0;
$prior_z = 0;
$prior_ats = 0;
#
#   ensure that every entity in the atom list has an entity number
#
$cur_model_count = 0;
$cur_model = '.';
for ($i = 1; $i < $atom_flag; ++$i) {
    while ($cur_model_count < $model_count &&	#???

      $i > $model_asn_high{$cur_model_count}) {	#???
	$cur_model_count++;
	$cur_model = $model{$cur_model_count};
    }
    if ($entity_seq_num{$entity_id{$i}} + 0 == 0) {
	$num_ents++;
	$entity_seq_num{$entity_id{$i}} = $num_ents;

	printf

	  "# *** WARNING *** At atom number %5s unidentified entity %3s found \n",

	  $atom_number{$i}, $entity_id{$i};

	$warning_list{++$warning_flag} =

	  sprintf("#=# ATOM_SITE: At atom number %5s unidentified entity %3s found \n",

	  $atom_number{$i}, $entity_id{$i});
    }
    #
    #  Assign an _entity_poly_seq.num value to each residue in the
    #  atom list
    #
    $atom_ent_seq_num{$i} = '.';
    if ($prior_res_name ne ($entity_id{$i} . $residue_seq_number{$i}) ||	#???

      $atom_or_het{$i} eq 'TER') {
	$ppoxt_found = $poxt_found;
	$poxt_found = $oxt_found;
	$oxt_found = 0;
	$ppprior_block = $pprior_block;
	$pprior_block = $prior_block;
	$prior_block = ($entity_id{$i} . '|' . $cur_model);
	$ppprior_res_name = $pprior_res_name;
	$pprior_res_name = $prior_res_name;
	$prior_res_name = ($entity_id{$i} . $residue_seq_number{$i});
	$pprior_x = $prior_x;
	$pprior_y = $prior_y;
	$pprior_z = $prior_z;
	$pprior_ats = $prior_ats;
	if ($al_res_ats > 0) {
	    $prior_x = $al_res_x / $al_res_ats;
	    $prior_y = $al_res_y / $al_res_ats;
	    $prior_z = $al_res_z / $al_res_ats;
	}
	$prior_ats = $al_res_ats;
	$at_bonded = '0';
	if ($prior_ats > 0 && $pprior_ats > 0 && $ppoxt_found == 0) {
	    $res_dist_sq = ($prior_x - $pprior_x) * ($prior_x - $pprior_x) +

	      ($prior_y - $pprior_y) * ($prior_y - $pprior_y) + ($prior_z -

	      $pprior_z) * ($prior_z - $pprior_z);
	    $res_max_sq = ($pprior_ats + $prior_ats) * ($pprior_ats +

	      $prior_ats) * 2;
	    if ($res_max_sq < 25) {
		$res_max_sq = 25;
	    }
	    $pnum = $pprior_ats + 1;
	    if ($prior_ats < $pnum) {	#???
		$pnum = $prior_ats + 1;
	    }
	    $res_min_sq = 2.5;
	    $at_bonded = '1';
	    if (($res_dist_sq < $res_min_sq) ||	#???

	      ($res_dist_sq > $res_max_sq)) {	#???
		$at_bonded = '0';
	    }
	}
	$al_res_x = 0;
	$al_res_y = 0;
	$al_res_z = 0;
	$al_res_ats = 0;
	{
	    $x_num = (@xxx = split(' ', $residue_name{$i}, 9999));
	    if ($x_num >= 1 && $xxx[$x_num] eq '' && ' ' eq ' ') {
		--$x_num;
	    }
	}
	$res_locs{$prior_res_name} = ($res_locs{$prior_res_name} . ' ' . $i);
	$cur_res = $xxx[1];
	if ($res_code{$cur_res} eq '') {
	    ++$num_res_name;
	    if ($num_res_name > length($charx)) {
		$num_res_name -= length($charx);
	    }
	    $res_code{$cur_res} = substr($charx, $num_res_name, 1);
	}
	if ($ppprior_res_name ne '' && $pprior_block ne '' &&

	  $pprior_block eq $ppprior_block) {	#???
	    $al_bond{$pprior_block} = ($al_bond{$pprior_block} . $at_bonded);
	}
	if ($atom_or_het{$i} ne 'TER') {
	    $al_seq{($entity_id{$i} . '|' . $cur_model)} =

	      ($al_seq{($entity_id{$i} . '|' . $cur_model)} .

	      $res_code{$cur_res});
	    $al_back_point{(length($al_seq{($entity_id{$i} . '|' .

	      $cur_model)}) . ($entity_id{$i} . '|' . $cur_model))} = $i;
	}
    }
    if ($atom_name{$i} eq ' OXT') {
	$oxt_found = 1;
    }
    if ($atom_or_het{$i} eq 'TER') {
	$oxt_found = 1;
	$prior_block = '';
	$pprior_block = '';
	$prior_res_name = '';
	$pprior_res_name = '';
	$prior_ats = 0;
	$al_res_ats = 0;
    }
    else {
	if ((substr($atom_type{$i}, 1, 2) ne ' H') &&

	  (substr($atom_type{$i}, 1, 2) ne ' D') &&

	  (substr($atom_type{$i}, 1, 2) ne ' T')) {
	    $al_res_x += $atom_x{$i};
	    $al_res_y += $atom_y{$i};
	    $al_res_z += $atom_z{$i};
	    $al_res_ats++;
	}
    }
    $al_seq_point{$i} = length($al_seq{($entity_id{$i} . '|' . $cur_model)});
}
if ($at_res_ats > 0) {
    $pprior_x = $prior_x;
    $pprior_y = $prior_y;
    $pprior_z = $prior_z;
    $pprior_ats = $prior_ats;
    $prior_x = $al_res_x / $al_res_ats;
    $prior_y = $al_res_y / $al_res_ats;
    $prior_z = $al_res_z / $al_res_ats;
    $prior_ats = $al_res_ats;
    $at_bonded = '0';
    if ($prior_ats > 0 && $pprior_ats > 0 && $poxt_found == 0) {
	$res_dist_sq = ($prior_x - $pprior_x) * ($prior_x - $pprior_x) +

	  ($prior_y - $pprior_y) * ($prior_y - $pprior_y) + ($prior_z -

	  $pprior_z) * ($prior_z - $pprior_z);
	$res_max_sq = ($pprior_ats + $prior_ats) * ($pprior_ats + $prior_ats)

	  * 2;
	if ($res_max_sq < 25) {
	    $res_max_sq = 25;
	}
	$pnum = $pprior_ats + 1;
	if ($prior_ats < $pnum) {	#???
	    $pnum = $prior_ats + 1;
	}
	$res_min_sq = 2.5;
	$at_bonded = '1';
	if (($res_dist_sq < $res_min_sq) || ($res_dist_sq > $res_max_sq)) {	#???	#???
	    $at_bonded = '0';
	}
    }
    if ($pprior_res_name ne '' && $prior_block ne '' &&

      $pprior_block eq $prior_block) {	#???
	$al_bond{$prior_block} = ($al_bond{$prior_block} . $at_bonded);
    }
}
foreach $iandm (keys %al_seq) {
    {
	$inum = (@xxiandm = split(/\|/, $iandm, 9999));
	if ($inum >= 1 && $xxiandm[$inum] eq '' && '|' eq ' ') {
	    --$inum;
	}
    }
    $i = $xxiandm[1];
    $mymodel = $xxiandm[2];
    $al_seq_c{$i} = '';
}
foreach $iandm (keys %al_seq) {
    {
	$inum = (@xxiandm = split(/\|/, $iandm, 9999));
	if ($inum >= 1 && $xxiandm[$inum] eq '' && '|' eq ' ') {
	    --$inum;
	}
    }
    $i = $xxiandm[1];
    $mymodel = $xxiandm[2];
    $nn_res = $entity_seq_num{$i};
    if ($al_seq_c{$i} ne '') {
	if ($al_seq{$iandm} ne $al_seq_c{$i}) {	#???
	    printf

	      "# *** WARNING *** sequence mismatch chain %3s, model %s\n", $i,

	      $mymodel;
	    $warning_list{++$warning_flag} =

	      sprintf("#=# ATOM_SITE: Sequence mismatch chain %3s, model %s\n",

	      int($xmat), $i, $mymodel);
	}
    }
    if ($nn_res <= $num_poly_ents && $al_seq_c{$i} eq '') {	#???
	$known_links = $al_bond{$iandm};
	{
	    $lsal = length($al_seq{$iandm});
	    $seqlim = $lsal / 2;
	    if ($seqlim > 6) {
		$seqlim = 6;
	    }
	    $lssr = length($seq_chain{$nn_res});
	    $isr = 1;
	    $nm = 0;
	    $missed = 1;
	    for ($isha = 1; $isha <= $lsal; ++$isha) {
		$seq_match{$isha} = 0;
	    }
	    for ($ipass = 0; (($ipass < 4) && ($missed != 0)); $ipass++) {
		$missed = 0;
		$isha = 1;
		$isr = 1;
		while ($isha <= $lsal) {	#???
		    $sm = 0;
		    $lm = $lsal - $isha + 1;
		    if ($seq_match{$isha} == 0) {
			while ($lm > 0 && $sm == 0) {
			    $cr = substr($al_seq{$iandm}, $isha, $lm);
			    $sm = index(substr($seq_chain{$nn_res}, $isr,

			      $lssr - $isr + 1), $cr);
			    $sf = 1;
			    if (!($known_links eq '') && !($sm == 0) &&

			      ($lm < $seqlim)) {	#???
				if ($isha < $lsal && $ipass > 0) {	#???
				    $kl = substr($known_links,

				      $isha + $lm - 1, 1);
				    if ($kl eq '1' && $ipass < 3 &&

				      $seq_match{$isha + $lm} > 0) {
					if (!($sm + $isr + $lm - 1 ==

					  $seq_match{$isha + $lm}) &&

					  ($ipass < 3)) {
					    $sm = 0;
					}
				    }
				    else {
					$kbase = $isha + $lm - 1;
					$kbreak = 0;
					while (($kbase < $lsal) &&	#???

					  !($seq_match{$kbase + 1} > 0)) {
					    $kbase++;
					    if (substr($known_links, $kbase,

					      1) eq '0') {
						$kbreak++;
					    }
					}
					if (($seq_match{$kbase + 1} > 0) &&

					  ((($sm + $isr - 1 + $kbase - $isha -

					  $lm > $seq_match{$kbase + 1} +

					  $kbreak * 3 - 1)) ||

					  ($sm + $isr - 1 + $kbase - $isha -

					  $lm < $seq_match{$kbase + 1} -

					  $kbreak * 3 - 2 - $seqlim * 2))) {
					    $sf = -1;
					}
				    }
				}
				if ($isha > 1 && !($sm == 0) && $sf == 1) {
				    $kl = substr($known_links, $isha, 1);
				    if ($kl eq '1' && $ipass < 3 &&

				      $seq_match{$isha - 1} > 0) {
					if (!($sm + $isr - 2 ==

					  $seq_match{$isha - 1}) &&

					  ($ipass < 3)) {
					    $sm = 0;
					}
				    }
				    else {
					$kbase = $isha;
					$kbreak = 0;
					while ((($kbase - 1) > 0) &&

					  !($seq_match{$kbase} > 0)) {
					    $kbase--;
					    if (substr($known_links, $kbase,

					      1) eq '0') {
						$kbreak++;
					    }
					}
					if (($seq_match{$kbase} > 0) &&

					  ((($sm + $isr - 1 + $kbase - $isha <

					  $seq_match{$kbase} - $kbreak * 3 -

					  1)) || $sm + $isr - 1 + $kbase -

					  $isha > $seq_match{$kbase} + $kbreak

					  * 3 + $seqlim * 2 - 1)) {
					    $sf = -1;
					}
				    }
				}
			    }

			    if ($sm > 0) {
				$kbase = $isha;
				while (($kbase > 1) &&

				  ($seq_match{$kbase - 1} <= 0) &&

				  ($sf > 0) &&

				  (substr($al_seq{$iandm}, $kbase - 1,

				  1) eq substr($seq_chain{$nn_res},	#???

				  $sm + $isr - 2 + $kbase - $isha, 1))) {
				    $kbase--;
				    $nm++;
				    $seq_match{$kbase} = $sm + $isr - 1 +

				      $kbase - $isha;
				}
				for ($ksa = 0; $ksa < length($cr); ++$ksa) {
				    if ($seq_match{$isha} == 0) {
					$seq_match{$isha} = ($sm + $ksa + $isr

					  - 1) * $sf;
					++$isha;
				    }
				    else {
					$cr = substr($cr, 1, $ksa);
				    }
				}
				$isr = $sm + length($cr) + $isr - 1;
				if (length($cr) > 3 || $isha == $lsal + 1) {
				    if ($sf > 0) {
					$nm += length($cr);
				    }
				}
				last;
			    }
			    else {
				if ($lm > 16) {
				    $lm = int($lm * .707);
				}
				else {
				    --$lm;
				    if ($ipass < 3 && $lm < $seqlim) {	#???
					if (($isha == 1) ||

					  ($isha > 1 &&

					  $seq_match{$isha - 1} == 0) &&

					  ($isha < $lsal - $lm + 1 &&

					  $seq_match{$isha + $lm} == 0)) {
					    $lm = 0;
					}
				    }
				}
			    }
			}
		    }
		    if ($sm == 0) {
			++$isha;
			$missed++;
		    }
		}
	    }
	    $qxmaxlen = $lsal;
	    if ($lssr > $lsal) {	#???
		$qxmaxlen = $lssr;
	    }
	    $ret_val = 100 * ($nm / $qxmaxlen);

	    $xmat = $ret_val;
	}

	if ($xmat < 90) {
	    printf "# *** WARNING *** only %5d%% homology to chain %3s\n",

	      int($xmat), $i;
	    $warning_list{++$warning_flag} =

	      sprintf("#=# ATOM_SITE: Only %5d%%%% homology to chain %3s\n",

	      int($xmat), $i);
	}
	$al_seq_c{$i} = $al_seq{$iandm};
	$first_seq_num{$i} = '';
	$last_seq_num{$i} = '';
	$misalign = 0;
	for ($j = 1; $j <= length($al_seq{$iandm}); ++$j) {
	    if ($seq_match{$j} > 0) {
		$k = $al_back_point{($j . $iandm)};
		$atom_ent_seq_num{$k} = $seq_match{$j} +

		  $seq_start_num{$nn_res} - 1;
		$ent_pol_seq_num{($i . $residue_seq_number{$k})} =

		  $atom_ent_seq_num{$k};
		if ($first_seq_num{$i} eq '') {
		    $first_seq_num{$i} = $atom_ent_seq_num{$k};
		}
		$last_seq_num{$i} = $atom_ent_seq_num{$k};
	    }
	}
	for ($j = 1; $j <= length($al_seq{$iandm}); ++$j) {
	    if ($seq_match{$j} <= 0) {
		$k = $al_back_point{($j . $iandm)};
		$prior_res_name = ($entity_id{$k} . $residue_seq_number{$k});
		$oepsn = $ent_pol_seq_num{$prior_res_name};
		if ($oepsn ne '' &&

		  ($oepsn = -$seq_match{$j} + $seq_start_num{$nn_res} - 1)) {
		    $atom_ent_seq_num{$k} = $oepsn;
		}
		else {
		    $misalign++;
		    printf

		      "# *** WARNING *** chain %3s residue %s %s not aligned to sequence\n",

		      $i, $residue_seq_number{$k}, $residue_name{$k};
		}
		$seq_match{$j} = 0;
	    }
	}
	if ($misalign > 0) {
	    $warning_list{++$warning_flag} =

	      sprintf("#=# ATOM_SITE: %d Unmatched residues in chain %3s\n",

	      $misalign, $i);
	}
    }
}
$prior_res_name = '';
for ($i = 1; $i < $atom_flag; ++$i) {
    if ($prior_res_name ne ($entity_id{$i} . $residue_seq_number{$i})) {	#???
	$prior_res_name = ($entity_id{$i} . $residue_seq_number{$i});
	$al_res_name{($chain_id{$i} . $residue_seq_number{$i})} =

	  $residue_name{$i};
	if ($atom_ent_seq_num{$i} eq '.') {
	    $oepsn = $ent_pol_seq_num{$prior_res_name};
	    if ($oepsn ne '') {
		$atom_ent_seq_num{$i} = $oepsn;
	    }
	}
    }
    if ($dense_list ne 'yes' && $atom_ent_seq_num{$i} eq '.') {
	$oepsn = $ent_pol_seq_num{$prior_res_name};
	if ($oepsn ne '') {
	    $atom_ent_seq_num{$i} = $oepsn;
	}
    }
}
#
#  Process DBREF, SEQADV, and MODRES information
#
$dbref_point{' '} = '';
$num_dbref = 0;
$num_align = 0;
$modres_chains = '';
$xdbref_flag = $dbref_flag;
for ($i = 1; $i <= $modres_flag; ++$i) {
    $mrc = $modres_chainID{$i};
    if (index($modres_chains, $mrc) == 0) {
	$modres_chains = ($modres_chain . $mrc);
	$dbref_database{++$dbref_flag} = '.';
	$dbref_chainID{$dbref_flag} = $mrc;
	$dbref_dbAccession{$dbref_flag} = '.';
	$dbref_dbIdCode{$dbref_flag} = '';
	$dbref_seqBegin{$dbref_flag} = '     ';
	$dbref_seqEnd{$dbref_flag} = '     ';
	$dbref_dbseqBeg{$dbref_flag} = '      ';
	$dbref_dbseqEnd{$dbref_flag} = '      ';
    }
}
for ($idb = 1; $idb <= $dbref_flag; ++$idb) {
    $dbref_align{$idb} = 'partial';
    ++$num_align;
    $dbn = $dbref_database{$idb};
    $cid = $dbref_chainID{$idb};
    $dbref_seq_dif{($cid . '|' . $dbn . '|' . $dbref_dbAccession{$idb})} =

      'no';
    $esn = $entity_seq_num{$cid};
    $dbref_esn{$idb} = $esn;
    $ssn = $seq_start_num{$esn};
    $fsn = $first_seq_num{$cid};
    $lsn = $last_seq_num{$cid};
    $dbfsn = '';
    if ($dbref_seqBegin{$idb} ne '     ') {
	$dbfsn = $ent_pol_seq_num{($cid . $dbref_seqBegin{$idb})};
	if ($dbfsn eq '') {	#???
	    $dbfsn = 1;
	}
    }
    $dblsn = '';
    if ($dbref_seqEnd{$idb} ne '     ') {
	$dblsn = $ent_pol_seq_num{($cid . $dbref_seqEnd{$idb})};
	if ($dblsn eq '') {	#???
	    $dblsn = length($seq_chain{$esn});
	}
    }
    $dbtfsn = '';
    $dbtlsn = '';
    if ($dbref_dbseqBeg{$idb} ne '      ' ||

      $dbref_dbseqEnd{$idb} ne '      ') {
	$dbtfsn = $dbref_dbseqBeg{$idb};
	$dbtlsn = $dbref_dbseqEnd{$idb};
	if ($dbref_database{$idb} eq 'PDB') {
	    $dbtfsn = $ent_pol_seq_num{($cid . $dbtfsn)};
	    if ($dbtfsn eq '') {	#???
		$dbtfsn = 1;
	    }
	    $dbtlsn = $ent_pol_seq_num{($cid . $dbtlsn)};
	    if ($dbtlsn eq '') {	#???
		$dbtlsn = length($seq_chain{$esn});
	    }
	}
	else {
	    $dbtfsn += 0;
	    $dbtlsn += 0;
	}
    }
    if ($dbfsn eq '' && $dblsn eq '' && $dbtfsn eq '' && $dbtlsn eq '') {	#???	#???	#???	#???
	$dbref_align{$idb} = '.';
	--$num_align;
	if ($modres_flag > 0 && $idb > $xdbref_flag) {	#???
	    $dbref_seq_dif{($cid . '|' . $dbn . '|' .

	      $dbref_dbAccession{$idb})} = 'yes';
	}
	else {
	    $dbref_esn{$idb} = '.';
	    $cid = $head_PDB_code;
	    $dbref_chainID{$idb} = $cid;
	    $dbref_seq_dif{($cid . '|' . $dbn . '|' .

	      $dbref_dbAccession{$idb})} = 'no';
	}
    }
    $dbid = ($dbref_dbAccession{$idb} . $dbref_dbIdCode{$idb});
    $dbref_tfsn{$idb} = $dbtfsn;
    $dbref_tlsn{$idb} = $dbtlsn;
    $dbref_fsn{$idb} = $dbfsn;
    $dbref_lsn{$idb} = $dblsn;
    $dbref_base_ref{$idb} = $idb;
    $prior_ref = ' ';
    $next_ref = $dbref_point{' '};
    $reloop = 'yes';
    while ($reloop eq 'yes') {
	if ($next_ref eq '') {
	    $dbref_point{$idb} = $next_ref;
	    $dbref_point{$prior_ref} = $idb;
	    ++$num_dbref;
	    $reloop = 'no';
	}
	else {
	    $xdbn = $dbref_database{$next_ref};
	    $xcid = $dbref_chainID{$next_ref};
	    $xdbid = ($dbref_dbAccession{$next_ref} .

	      $dbref_dbIdCode{$next_ref});
	    if ($cid lt $xcid || ($cid eq $xcid && $dbn lt $xdbn) ||	#???	#???	#???

	      ($cid eq $xcid && $dbn eq $xdbn && $dbid lt $xdbid)) {	#???	#???	#???
		$dbref_point{$idb} = $next_ref;
		$dbref_point{$prior_ref} = $idb;
		++$num_dbref;
		$reloop = 'no';
	    }
	    else {
		if ($cid eq $xcid && $dbn eq $xdbn && $dbid eq $xdbid) {	#???	#???	#???
		    $dbref_base_ref{$idb} = $next_ref;
		    $reloop = 'no';
		}
		else {
		    $prior_ref = $next_ref;
		    $next_ref = $dbref_point{$prior_ref};
		}
	    }
	}
    }
}
#
#  Merge in information about point differences
#
for ($isq = 1; $isq <= $seqadv_flag; ++$isq) {
    $dbn = $seqadv_database{$isq};
    $cid = $seqadv_chainID{$isq};
    $dbid = $seqadv_dbAccession{$isq};
    if ($dbref_seq_dif{($cid . '|' . $dbn . '|' . $dbid)} eq '') {
	$warning_list{++$warning_flag} =

	  ('#=# STRUCT_REF: Unable to locate DBREF for ' . $cid . ' ' . $dbn .

	  ' ' . $dbid . "\n");
    }
    $dbref_seq_dif{($cid . '|' . $dbn . '|' . $dbid)} = 'yes';
}

#
#  If there are multiple models, verify the alignment of atoms
#
if ($model_flags eq 'yes') {
    #
    # prepare pointers for residue list
    #
    $last_alres = 0;
    $first_alres = 0;
    for ($imodel = 1; $imodel <= $model_count; $imodel++) {
	$jmodel = $imodel + 1;
	if ($jmodel == $model_count + 1) {
	    $jmodel = 1;
	}
	$iaflo = $model_asn_low{$imodel};
	$iafhi = $model_asn_high{$imodel};
	$jaflo = $model_asn_low{$jmodel};
	$jafhi = $model_asn_high{$jmodel};
	if ($iafhi - $iaflo != $jafhi - $jaflo) {
	    $warning_list{++$warning_flag} = ('#=# ATOM_SITE: Model ' .

	      $model{$imodel} . " inconsistent number of atoms \n");
	}
	$imatch = 'true';
	$kerror = 0;
	$iasn = $iaflo;
	$jasn = $jaflo;
	while ($iasn <= $iafhi) {	#???
	    if ($jasn > $jafhi || $atom_name{$iasn} ne $atom_name{$jasn} ||	#???	#???

	      $residue_name{$iasn} ne $residue_name{$jasn} ||	#???

	      $atom_alt_location{$iasn} ne $atom_alt_location{$jasn} ||	#???

	      $chain_id{$iasn} ne $chain_id{$jasn} ||	#???

	      $residue_seq_number{$iasn} ne $residue_seq_number{$jasn}) {	#???
		$imatch = 'false';
		$kerror++;
	    }
	    else {
		if ($root_model_rep{$iasn} eq '' && $imodel == 1) {
		    $root_model_rep{$iasn} = $iasn;
		}
		if ($root_model_rep{$jasn} eq '') {
		    $root_model_rep{$jasn} = $root_model_rep{$iasn};
		}
	    }
	    $alres = ($residue_name{$iasn} . '|' . $chain_id{$iasn} . '|' .

	      $residue_seq_number{$iasn});
	    $kpoint = $reslist{$alres};
	    if ($kpoint eq '') {
		$next_alres{$last_alres} = $iasn;
		$last_alres = $iasn;
		$reslist{$alres} = $iasn;
		$next_atom{$iasn} = '';
		if ($root_model_rep{$iasn} eq '') {
		    $root_model_rep{$iasn} = $iasn;
		}
	    }
	    else {
		$alatom = ($atom_name{$iasn} . '|' .

		  $atom_alt_location{$iasn});
		while ($kpoint ne '') {
		    $blatom = ($atom_name{$kpoint} . '|' .

		      $atom_alt_location{$kpoint});
		    if ($alatom ne $blatom) {	#???
			$lpoint = $kpoint;
			$kpoint = $next_atom{$kpoint};
			if ($kpoint eq '') {
			    $next_atom{$lpoint} = $iasn;
			    $next_atom{$iasn} = '';
			    if ($imodel - 1 != 0) {
				$warning_list{++$warning_flag} =

				  ('#=# ATOM_SITE: Model ' . $model{$imodel} .

				  ' added new atom ' . $iasn . ' ' .

				  $atom_number{$iasn} . "\n");
			    }
			    if ($root_model_rep{$iasn} eq '') {
				$root_model_rep{$iasn} = $iasn;
			    }
			}
		    }
		    else {
			if ($imodel == 1) {
			    $warning_list{++$warning_flag} =

			      ('#=# ATOM_SITE: Atom ' . $iasn .

			      ' duplicates atom ' . $kpoint . "\n");
			}
			$root_model_rep{$iasn} = $kpoint;
			$kpoint = '';
		    }
		}
	    }
	    $iasn++;
	    $jasn++;
	}
	if ($imatch ne 'true') {
	    $warning_list{++$warning_flag} = ('#=# ATOM_SITE: Model ' .

	      $model{$imodel} . ' fails to match ' . $kerror .

	      " atoms to next model \n");
	}
    }
    $respoint = $next_alres{$first_alres};
    $id_in_model = 0;
    $max_atom_number = '';
    while ($respoint + 0 != 0) {
	$iasn = $respoint;
	$respoint = $next_alres{$respoint};
	while ($iasn != 0) {
	    $atom_id_by_type{($atom_name{$iasn} . '|' . $residue_name{$iasn} .

	      '|' . $atom_alt_location{$iasn} . '|' . $chain_id{$iasn} . '|' .

	      $residue_seq_number{$iasn})} = ++$id_in_model;
	    if ($root_model_rep{$iasn} eq '') {
		$root_model_rep{$iasn} = $iasn;
	    }
	    $atom_id_in_model{$root_model_rep{$iasn}} = $id_in_model;
	    $atom_point_id_in_model{$id_in_model} = $iasn;
	    if ($max_atom_number eq '') {	#???
		$max_atom_number = $atom_number{$iasn};
	    }
	    else {
		if ($max_atom_number - $atom_number{$iasn} < 0) {
		    $max_atom_number = $atom_number{$iasn};
		}
	    }
	    $iasn = $next_atom{$iasn} + 0;
	}
    }
    $max_asn = $atom_point{$max_atom_number + 0};
}

printf (("\nloop_\n"));
printf (("_entity.id\n"));
printf (("_entity.type\n"));
printf (("_entity.details\n"));
for ($i = 1; $i <= $num_poly_ents; ++$i) {
    if ($ent_poly_num{$entities{$i}} eq $entity_seq_num{$entities{$i}}) {	#???
	printf " %5d     polymer\n; %s\n;\n", $entity_seq_num{$entities{$i}},

	  $entities_list{$i};
    }
}
for ($i = $num_poly_ents + 1; $i < $num_ents; ++$i) {
    if ($entities{$i} ne 'HOH' && $entities{$i} ne 'DOD') {
	printf " %5d     non-polymer 'het group %s'\n",

	  $entity_seq_num{$entities{$i}}, $entities{$i};
    }
    else {
	printf " %5d     water       '%s' \n", $entity_seq_num{$entities{$i}},

	  $entities{$i};
    }
}

if ($hetnam_flag + 0 > 1) {
    printf (("\nloop_\n"));
    printf (("_entity_name_com.entity_id\n"));
    printf (("_entity_name_com.name\n\n"));
    foreach $hnames (keys %het_site_name) {
	$num_hnames = $entity_seq_num{$hnames};
	if ($num_hnames eq '') {
	    $warning_list{++$warning_flag} =

	      ('#=# ENTITY_NAME:  Unable to locate entity.id for ' . $hnames .

	      "\n");
	    $num_hnames = '.';
	}
	if ($num_hnames ne '.') {
	    printf "%5d\n;%s\n;\n", $num_hnames, $het_site_name{$hnames};
	}
	else {
	    printf "  %s  \n;%s\n;\n", $num_hnames, $het_site_name{$hnames};
	}
	$hsyn = $het_site_syn{$hnames};
	if ($hsyn ne '') {
	    {
		$num_hsyn = (@hsyns = split(/;/, $hsyn, 9999));
		if ($num_hsyn >= 1 && $hsyns[$num_hsyn] eq '' && ';' eq ' ') {
		    --$num_hsyn;
		}
	    }
	    for ($i = 1; $i <= $num_hsyn; ++$i) {
		$xhsyn = $hsyns[$i];
		{
		    $num_xhsyn = (@xhsyns = split(' ', $xhsyn, 9999));
		    if ($num_xhsyn >= 1 && $xhsyns[$num_xhsyn] eq '' &&

		      ' ' eq ' ') {
			--$num_xhsyn;
		    }
		}
		$yhsyn = '';
		$ll = 0;
		for ($ii = 1; $ii <= $num_xhsyn; ++$ii) {
		    if (length($xhsyns[$ii]) + $ll > 77) {
			if ($ll > 0) {
			    $yhsyn = ($yhsyn . "\n  ");
			}
			else {
			    $yhsyn = ($yhsyn . substr($xhsyns[$ii], 1,

			      77) . "\n  ");
			    $xhsyns[$ii] = substr($xhsyns[$ii], 78,

			      length($xhsyns[$ii] - 77));
			}
			$ll = 0;
			--$ii;
		    }
		    else {
			if ($ll > 0) {
			    $yhsyn = ($yhsyn . ' ' . $xhsyns[$ii]);
			    $ll = $ll + 1 + length($xhsyns[$ii]);
			}
			else {
			    $yhsyn = ($yhsyn . $xhsyns[$ii]);
			    $ll = length($xhsyns[$ii]);
			}
		    }
		}
		if (length($yhsyn) > 0) {
		    if ($num_hnames ne '.') {
			printf "%5d\n; %s\n;\n", $num_hnames, $yhsyn;
		    }
		    else {
			printf "  %s  \n; %s\n;\n", $num_hnames, $yhsyn;
		    }
		}
	    }
	}
    }
}

printf (("\nloop_\n"));
printf (("_struct_asym.entity_id\n"));
printf (("_struct_asym.id\n"));
for ($i = 1; $i < $num_ents; ++$i) {
    printf " %5d %s\n", $entity_seq_num{$entities{$i}}, $entities{$i};
}
if ($verbose eq 'yes') {
    printf (("\nloop_\n"));
    printf (("_entity_keywords.entity_id\n"));
    printf (("_entity_keywords.text\n"));
    for ($i = 1; $i <= $num_poly_ents; ++$i) {
	printf " %3s     ? \n", $entity_seq_num{$entities{$i}};
    }
    for ($i = $num_poly_ents + 1; $i < $num_ents; ++$i) {
	printf " %3s     ? \n", $entity_seq_num{$entities{$i}};
    }
}
#
# print out DBREF information
#
if ($num_dbref + 0 > 0) {
    printf (("\nloop_\n"));
    printf (("_struct_ref.id \n"));
    printf (("_struct_ref.entity_id \n"));
    printf (("_struct_ref.biol_id \n"));
    printf (("_struct_ref.db_name \n"));
    printf (("_struct_ref.db_code \n"));
    printf (("_struct_ref.seq_align \n"));
    printf (("_struct_ref.seq_dif \n"));
    printf (("_struct_ref.details \n"));
    $xnumdbref = 0;
    $idb = ' ';
    for ($i = 1; $i <= $num_dbref; ++$i) {
	$idb = $dbref_point{$idb};
	++$xnumdbref;
	$dbref_id{($dbref_chainID{$idb} . '|' . $dbref_database{$idb} . '|' .

	  $dbref_dbAccession{$idb})} = $xnumdbref;
	$esn = $dbref_esn{$idb};
	$dbac = $dbref_dbAccession{$idb};
	if ($dbref_dbIdCode{$idb} ne '') {
	    $dbac = ($dbref_dbAccession{$idb} . $dbref_dbIdCode{$idb});
	}
	if ($dbac ne '.') {
	    $dbac = ("'" . $dbac . "'");
	}
	$comment = '.';
	printf " %5d %4s %4s %8s  %s  %7s  %3s %s\n", $xnumdbref, $esn,

	  $dbref_chainID{$idb}, $dbref_database{$idb}, $dbac,

	  $dbref_align{$idb},

	  $dbref_seq_dif{($dbref_chainID{$idb} . '|' . $dbref_database{$idb} .

	  '|' . $dbref_dbAccession{$idb})}, $comment;
    }
    printf (("\nloop_\n"));
    printf (("_struct_ref_seq.align_id \n"));
    printf (("_struct_ref_seq.ref_id \n"));
    printf (("_struct_ref_seq.seq_align_beg \n"));
    printf (("_struct_ref_seq.seq_align_end\n"));
    printf (("_struct_ref_seq.db_align_beg \n"));
    printf (("_struct_ref_seq.db_align_end \n"));
    printf (("_struct_ref_seq.details \n"));
    #
    #  Print out alignment information for each chain fragment
    #
    for ($idb = 1; $idb <= $dbref_flag; ++$idb) {
	$dbn = $dbref_database{$idb};
	$cid = $dbref_chainID{$idb};
	$xnumdbref = $dbref_id{($dbref_chainID{$idb} . '|' .

	  $dbref_database{$idb} . '|' . $dbref_dbAccession{$idb})};
	$dbref_align_id{($dbref_chainID{$idb} . '|' . $dbref_database{$idb} .

	  '|' . $dbref_dbAccession{$idb})} = $idb;
	$dbfsn = '   .   ';
	$dblsn = '   .   ';
	$dbtfsn = '   .   ';
	$dbtlsn = '   .   ';
	if ($dbref_fsn{$idb} ne '') {
	    $dbfsn = $dbref_fsn{$idb};
	}
	if ($dbref_lsn{$idb} ne '') {
	    $dblsn = $dbref_lsn{$idb};
	}
	if ($dbref_tfsn{$idb} ne '') {
	    $dbtfsn = $dbref_tfsn{$idb};
	}
	if ($dbref_tlsn{$idb} ne '') {
	    $dbtlsn = $dbref_tlsn{$idb};
	}
	$dbid = ($dbref_dbAccession{$idb} . $dbref_dbIdCode{$idb});
	if ($dbref_align{$idb} eq 'partial' ||

	  ($dbref_align{$idb} eq '.' && $idb > $xdbref_flag)) {	#???
	    printf " %5d %5d %7s %7s %7s %7s . \n", $idb, $xnumdbref, $dbfsn,

	      $dblsn, $dbtfsn, $dbtlsn;
	}
    }
    if ($seqadv_flag > 0 || $modres_flag > 0) {
	printf (("\nloop_\n"));
	printf (("_struct_ref_seq_dif.align_id \n"));
	printf (("_struct_ref_seq_dif.seq_num \n"));
	printf (("_struct_ref_seq_dif.mon_id \n"));
	printf (("_struct_ref_seq_dif.db_seq_num \n"));
	printf (("_struct_ref_seq_dif.db_mon_id \n"));
	printf (("_struct_ref_seq_dif.details \n"));
	for ($isq = 1; $isq <= $seqadv_flag; ++$isq) {
	    $dbn = $seqadv_database{$isq};
	    $cid = $seqadv_chainID{$isq};
	    $esn = $entity_seq_num{$cid};
	    $seq_num = $seqadv_seq{$isq};
	    if ($seq_num eq '') {
		$seq_num = '  .  ';
	    }
	    $resName = $seqadv_resName{$isq};
	    $dbid = $seqadv_dbAccession{$isq};
	    $db_seq_num = $seqadv_dbSeq{$isq};
	    if ($db_seq_num eq '') {
		$db_seq_num = '  .  ';
	    }
	    $db_resName = $seqadv_dbRes{$isq};
	    $comment = $seqadv_conflict{$isq};
	    printf (("%5d %5s %3s %5s %3s\n" . "         %s\n"),

	      $dbref_align_id{($cid . '|' . $dbn . '|' .

	      $seqadv_dbAccession{$isq})}, $seq_num, $resName, $db_seq_num,

	      $db_resName, $comment);
	}
	for ($isq = 1; $isq <= $modres_flag; ++$isq) {
	    $cid = $modres_chainID{$isq};
	    $esn = $entity_seq_num{$cid};
	    $seq_num = $modres_seq{$isq};
	    if ($seq_num eq '' || $seq_num eq '     ') {
		$seq_num = '  .  ';
	    }
	    $resName = $modres_resName{$isq};
	    $db_resName = $modres_dbRes{$isq};
	    $comment = $modres_conflict{$isq};
	    printf (("%5d %5s %3s %5s %3s\n" . "         %s\n"),

	      $dbref_align_id{($cid . '|.|.')}, $seq_num, $resName, '.',

	      $db_resName, $comment);
	}
    }
}
if ($verbose eq 'yes') {
    if ($num_dbref + 0 == 0) {
	printf (("\nloop_\n"));
	printf (("_struct_ref.entity_id \n"));
	printf (("_struct_ref.db_name \n"));
	printf (("_struct_ref.db_code \n"));
	printf (("_struct_ref.details \n"));
	for ($i = 1; $i <= $num_poly_ents; ++$i) {
	    if ($ent_poly_num{$entities{$i}} eq	#???

	      $entity_seq_num{$entities{$i}}) {
		printf " %3s     ?  ?  ? \n", $entity_seq_num{$entities{$i}};
	    }
	}
	for ($i = $num_poly_ents + 1; $i < $num_ents; ++$i) {
	    printf " %3s     ?  ?  ? \n", $entity_seq_num{$entities{$i}};
	}
    }
}
if ($verbose eq 'yes') {
    #
    # Information about entities which are polymers
    #
    printf (("\nloop_ \n"));
    printf (("_entity_poly.entity_id\n"));
    printf (("_entity_poly.type\n"));
    printf (("_entity_poly.nstd_chirality\n"));
    printf (("_entity_poly.nstd_linkage\n"));
    printf (("_entity_poly.nstd_monomer\n"));
    printf (("_entity_poly.type_details\n"));

    for ($i = 1; $i <= $num_poly_ents; ++$i) {
	printf "  %1s  ?   ?   ?   ?   ?\n", $entity_seq_num{$entities{$i}};
    }
}

#   Non-standard monomer entities described by FORMUL records
if ($formul_flag > 1 || $atom_res_flag > 0 || $res_res_flag > 0) {
    printf (("\n\n\n"));
    printf (("####################\n"));
    printf (("#                  #\n"));
    printf (("# CHEM_COMP        #\n"));
    printf (("#                  #\n"));
    printf (("####################\n\n\n"));

    printf (("loop_\n"));
    printf (("_chem_comp.id\n"));
    printf (("_chem_comp.mon_nstd_flag\n"));
    printf (("_chem_comp.formula\n"));
    printf (("_chem_comp.name\n"));
    if ($hetnam_flag + 0 > 1) {
	printf (("_chem_comp.mon_nstd_details\n"));
    }
    for ($i = 1; $i < $formul_flag; ++$i) {
	if ($formul_het_cont_flag{$i} eq '  ') {
	    $x_cur_res = substr(('   ' . $formul_het_site_symbol{$i}),

	      1 + length($formul_het_site_symbol{$i}), 3);
	    $res_count{$x_cur_res} = 0;
	    $mon_ns = 'no ';
	    $het_details = $het_site_text{$formul_het_site_symbol{$i}};
	    if ($het_details eq '') {
		$het_details = ' ';
	    }
	    if ($hetnam_flag + 0 > 1) {
		$het_name = $het_site_name{$formul_het_site_symbol{$i}};
		if ($het_name eq '') {
		    $het_name = ' ';
		}
		printf "%s %s%s\n;\n; %s\n;\n; %s\n;\n", 
		$formul_het_site_symbol{$i}, $mon_ns,

		  $het_formula{$formul_het_site_symbol{$i}}, $het_name,

		  $het_details;
	    }
	    else {
		printf "%s %s%s\n;\n; %s\n;\n", 
		$formul_het_site_symbol{$i}, $mon_ns,

		  $het_formula{$formul_het_site_symbol{$i}}, $het_details;
	    }
	}
    }
    # end for (i=1; i < formul_flag; ++i)
    for ($i = 1; $i <= $num_AARESES + $num_NARESES + $num_OTRESES; ++$i) {
	$x_res_code = $res_id{$i};
	if ($res_count{$x_res_code} + 0 > 0) {
	    $mon_formul = '.';
	    if ($res_formul{$x_res_code} ne '.') {
		$mon_formul = ("'" . $res_formul{$x_res_code} . "'");
	    }
	    $mon_name = '.';
	    if ($res_name{$x_res_code} ne '.') {
		$mon_name = ("\"" . $res_name{$x_res_code} . "\"");
	    }
	    $mon_ns = 'yes';
	    $mon_det = ' ';
	    if ($hetnam_flag + 0 > 1) {
		$mon_det = '.';
	    }
	    if ($i > $num_AARESES + $num_NARESES) {
		$mon_ns = 'no';
	    }
	    printf "%3s %3s %-25s %-25s %s\n", 
	    $x_res_code, $mon_ns, $mon_formul, $mon_name, $mon_det;
	}
    }
}

printf (("\n\n\n"));
printf (("######################\n"));
printf (("#                    #\n"));
printf (("# ATOM_SITES         #\n"));
printf (("#                    #\n"));
printf (("######################\n\n\n"));

# PRINT ORIGX INFO

for ($i = 1; $i <= $origx_flag; ++$i) {
    printf (($om_save{$i}));
}

# PRINT SCALE INFO

for ($i = 1; $i <= $scale_flag; ++$i) {
    printf (($sc_save{$i}));
}

# PRINT MTRIX INFO

for ($i = 1; $i <= $mtrix_flag; ++$i) {
    printf (($mat_save{$i}));
}

# PRINT TVECT INFO

for ($i = 1; $i <= $tvect_flag; ++$i) {
    printf (($tv_save{$i}));
}

#   print alternate sites info
if ($at_alt) {
    printf (("\n\n\n"));
    printf (("####################\n"));
    printf (("#                  #\n"));
    printf (("# ATOM_SITES_ALT   #\n"));
    printf (("#                  #\n"));
    printf (("####################\n\n\n"));

    printf (("\nloop_\n"));
    printf (("_atom_sites_alt.id\n"));
    printf (("_atom_sites_alt.details\n"));
    foreach $i (keys %atom_alt_list) {
	printf " %3s    ? \n", $i;
    }
    if ($verbose eq 'yes') {
	printf (("\n\n\n"));
	printf (("######################\n"));
	printf (("#                    #\n"));
	printf (("# ATOM_SITES_ALT_ENS #\n"));
	printf (("#                    #\n"));
	printf (("######################\n\n\n"));

	printf (("\nloop_\n"));
	printf (("_atom_sites_alt_ens.id\n"));
	printf (("_atom_sites_alt_ens.details\n"));
	printf ((" 'Ensemble 1'   ? \n"));

	printf (("\n\n\n"));
	printf (("######################\n"));
	printf (("#                    #\n"));
	printf (("# ATOM_SITES_ALT_GEN #\n"));
	printf (("#                    #\n"));
	printf (("######################\n\n\n"));

	printf (("\nloop_\n"));
	printf (("_atom_sites_alt_gen.ens_id\n"));
	printf (("_atom_sites_alt_gen.alt_id\n"));
	foreach $i (keys %atom_alt_list) {
	    printf " 'Ensemble 1'        %3s    \n", $i;
	}
    }
    # end if (verbose == "yes")

    ;
}
# end if (at_alt)

if ($foot_flag > 0) {
    printf (("\n\n\n"));
    printf (("######################\n"));
    printf (("#                    #\n"));
    printf (("# ATOM_SITES_FOOTNOTE#\n"));
    printf (("#                    #\n"));
    printf (("######################\n\n\n"));
    for ($i = 1; $i <= $foot_flag; ++$i) {
	printf (($ft_save{$i}));
    }
    printf ((";\n"));
}

#   print ATOM/HETATM info
if ($atom_flag_1 eq '1' && $atom_flag + 0 > 1) {	#???
    printf (("\n\n\n"));
    printf (("####################\n"));
    printf (("#                  #\n"));
    printf (("# ATOM_SITE        #\n"));
    printf (("#                  #\n"));
    printf (("####################\n\n\n"));

    printf (("\nloop_\n"));
    printf (("_atom_site.label_seq_id\n"));
    if ($model_flags eq 'yes') {
	printf ((('_atom_site.' . $pdb2cif_prefix . "label_model_id\n")));
	printf ((('_atom_site.' . $pdb2cif_prefix . "id_in_model\n")));
	printf ((('_atom_site.' . $pdb2cif_prefix . "auth_id_in_model\n")));
	if ($verbose eq 'yes') {
	    printf ((('_atom_site.' . $pdb2cif_prefix . "auth_model_id\n")));
	}
    }
    if ($compliance_level >= 2.0) {
	printf (("_atom_site.auth_asym_id\n"));
    }
    printf (("_atom_site.group_PDB\n"));
    printf (("_atom_site.type_symbol \n"));
    printf (("_atom_site.label_atom_id \n"));
    printf (("_atom_site.label_comp_id \n"));
    printf (("_atom_site.label_asym_id \n"));
    printf (("_atom_site.auth_seq_id \n"));
    printf (("_atom_site.label_alt_id \n"));
    printf (("_atom_site.cartn_x \n"));
    printf (("_atom_site.cartn_y \n"));
    printf (("_atom_site.cartn_z \n"));
    printf (("_atom_site.occupancy\n"));
    printf (("_atom_site.B_iso_or_equiv \n"));
    printf (("_atom_site.footnote_id\n"));
    printf (("_atom_site.label_entity_id\n"));
    printf (("_atom_site.id\n"));

    if ($sigatm_flag > 0) {
	printf (("_atom_site.cartn_x_esd \n"));
	printf (("_atom_site.cartn_y_esd \n"));
	printf (("_atom_site.cartn_z_esd \n"));
	printf (("_atom_site.occupancy_esd\n"));
	printf (("_atom_site.B_iso_or_equiv_esd \n"));
    }

    if ($aniso_flag > 0) {
	printf (("_atom_site.aniso_U[1][1]\n"));
	printf (("_atom_site.aniso_U[1][2]\n"));
	printf (("_atom_site.aniso_U[1][3]\n"));
	printf (("_atom_site.aniso_U[2][2]\n"));
	printf (("_atom_site.aniso_U[2][3]\n"));
	printf (("_atom_site.aniso_U[3][3]\n"));
    }
    # end if (aniso_flag > 0)

    if ($siguij_flag > 0) {
	printf (("_atom_site.aniso_U[1][1]_esd\n"));
	printf (("_atom_site.aniso_U[1][2]_esd\n"));
	printf (("_atom_site.aniso_U[1][3]_esd\n"));
	printf (("_atom_site.aniso_U[2][2]_esd\n"));
	printf (("_atom_site.aniso_U[2][3]_esd\n"));
	printf (("_atom_site.aniso_U[3][3]_esd\n"));
    }
    # end if (siguij_flag > 0)

    $atom_flag_1 = 2;
}
# end if (atom_flag_1 == "1")
$b_warn = 0;
$occ_warn = 0;
$cur_model = '.';
$cur_model_count = 0;
for ($i = 1; $i < $atom_flag; ++$i) {
    while ($cur_model_count < $model_count &&	#???

      $i > $model_asn_high{$cur_model_count}) {	#???
	$cur_model_count++;
	$cur_model = $model{$cur_model_count};
    }
    $es = ' ';
    $xes = '';
    $ts = ' ';
    if (substr($atom_type{$i}, 3, 2) ne '') {
	$ts = "\n          ";
    }
    if ($atom_ent_seq_num{$i} ne '.' || $model_flags eq 'yes' ||

      $compliance_level >= 2.0) {
	$es = "\n  ";
	$xes = '';
    }
    $xs = '';
    if (substr($atom_x{$i}, 1, 1) ne ' ') {
	$xs = "\n                             ";
    }
    $ys = '';
    if (substr($atom_y{$i}, 1, 1) ne ' ') {
	$ys = "\n                                     ";
    }
    $zs = '';
    if (substr($atom_z{$i}, 1, 1) ne ' ') {
	$zs = "\n                                             ";
    }
    $os = '';
    if (substr($atom_occ{$i}, 1, 1) ne ' ') {
	$os = "\n                                                     ";
    }
    $bs = '';
    if (substr($B_or_U{$i}, 1, 1) ne ' ') {
	$bs = "\n                                                           ";
    }
    $aco = '';
    if ($atom_or_het{$i} eq 'TER' && $print_ter eq 'comment') {
	$aco = '#';
	$ts = ' ';
	if (substr($atom_type{$i}, 3, 2) ne '') {
	    $ts = "\n#         ";
	}
	if ($compliance_level >= 2.0) {
	    $es = ' ';
	    $xes = "\n#    ";
	}
	else {
	    $es = "\n# ";
	    $xes = '';
	}
	if ($xs ne '') {
	    $xs = "\n#                            ";
	}
	if ($ys ne '') {
	    $ys = "\n#                                    ";
	}
	if ($zs ne '') {
	    $zs = "\n#                                            ";
	}
	if ($os ne '') {
	    $os = "\n#                                                    ";
	}
	if ($bs ne '') {
	    $bs =

	      "\n#                                                          ";
	}
    }
    if ($cur_model eq '') {
	$cur_model = '.';
    }
    if ($compliance_level >= 2.0) {
	$ts = '   ';
	if (substr($atom_type{$i}, 3, 2) ne '') {
	    $ts = '  ';
	}
	$asi = ("'" . $atom_seg_id{$i} . "'");
	if ($atom_seg_id{$i} eq '    ') {
	    $asi = '.';
	}
	delete $atom_seq_id{$i};
	$es = (' ' . $asi . "\n");
	$xes = '';
	if ($atom_or_het{$i} eq 'TER' && $print_ter eq 'comment') {
	    $es = (' ' . $asi . ' ');
	    $xes = "\n#    ";
	}
    }
    $anum = $atom_number{$i};
    if ($model_flags eq 'yes') {
	$id_in_model = 0;
	if ($root_model_rep{$i} ne '') {
	    $id_in_model = $atom_id_in_model{$root_model_rep{$i}};
	    $auth_id_in_model = $atom_number{$root_model_rep{$i}};
	}
	if ($id_in_model + 0 == 0) {
	    $id_in_model = $atom_id_by_type{($atom_name{$i} . '|' .

	      $residue_name{$i} . '|' . $atom_alt_location{$i} . '|' .

	      $chain_id{$i} . '|' . $residue_seq_number{$i})};
	    $auth_id_in_model =

	      $atom_number{$atom_point_id_in_model{$id_in_model}};
	}
	if ($id_in_model + 0 == 0) {
	    $id_in_model = '.';
	    $auth_id_in_model = '.';
	}
	if ($verbose ne 'yes') {
	    $es = (' ' . $cur_model . ' ' . $id_in_model . ' ' .

	      $auth_id_in_model . $es);
	}
	else {
	    $es = (' ' . $cur_model . ' ' . $id_in_model . ' ' .

	      $auth_id_in_model . ' .' . $es);
	}
	if ($atom_unique_ids eq 'no') {
	    $anum = ($cur_model_count - 1) * $max_atom_number +

	      $auth_id_in_model;
	}
    }
    $ans = ' ';
    if ($anum - 999999 > 0 && $anum - 9999999 <= 0) {
	$ans =

	  "\n                                                                         ";
	if ($atom_or_het{$i} eq 'TER' && $print_ter eq 'comment') {
	    $ans =

	      "\n#                                                                        ";
	}
    }
    else {
	if ($anum - 9999999 > 0 && $anum - 99999999 <= 0) {
	    $ans =

	      "\n                                                                        ";
	    if ($atom_or_het{$i} eq 'TER' && $print_ter eq 'comment') {
		$ans =

		  "\n#                                                                       ";
	    }
	}
	else {
	    if ($anum - 99999999 > 0) {
		$ans =

		  "\n                                                                       ";
		if ($atom_or_het{$i} eq 'TER' && $print_ter eq 'comment') {
		    $ans =

		      "\n#                                                                      ";
		}
	    }
	}
    }
    if ($atom_or_het{$i} ne 'TER' || $print_ter ne 'no') {
	if ($compliance_level >= 2.0) {
	    printf

	      "%s%s%-4s %s%s%s%4s %4s %1s %5s %1s%s%8s%s%8s%s%8s%s%6s%s%6s %s %3d%s%6d\n",

	      ($aco . $atom_ent_seq_num{$i}), $es, substr($atom_or_het{$i}, 1,

	      4), $xes, $atom_type{$i}, $ts, $atom_name{$i},

	      $residue_name{$i}, $chain_id{$i}, $residue_seq_number{$i},

	      $atom_alt_location{$i}, $xs, $atom_x{$i}, $ys, $atom_y{$i}, $zs,

	      $atom_z{$i}, $os, $atom_occ{$i}, $bs, $B_or_U{$i},

	      $footnote_number{$i}, $entity_seq_num{$entity_id{$i}}, $ans,

	      $anum;
	}
	else {
	    printf

	      "%s%s%-4s %s%s%4s %4s %1s %5s %1s%s%8s%s%8s%s%8s%s%6s%s%6s %s %3d%s%6d\n",

	      ($aco . $atom_ent_seq_num{$i}), $es, substr($atom_or_het{$i}, 1,

	      4), $atom_type{$i}, $ts, $atom_name{$i}, $residue_name{$i},

	      $chain_id{$i}, $residue_seq_number{$i}, $atom_alt_location{$i},

	      $xs, $atom_x{$i}, $ys, $atom_y{$i}, $zs, $atom_z{$i}, $os,

	      $atom_occ{$i}, $bs, $B_or_U{$i}, $footnote_number{$i},

	      $entity_seq_num{$entity_id{$i}}, $ans, $anum;
	}
	if ($B_or_U{$i} ne '   .  ') {
	    if ($B_or_U{$i} + 0 > 0) {
		++$b_warn;
	    }
	}
	if ($atom_occ{$i} ne '   .  ') {
	    if ($atom_occ{$i} + 0 < 0 || $atom_occ{$i} + 0 > 1) {
		++$occ_warn;
	    }
	    if ($nmr_flag > 0) {
		if ($atom_occ{$i} + 0 != 0 && $atom_occ{$i} + 0 != 1) {
		    ++$occ_warn;
		}
	    }
	}
	delete $atom_x{$i};
	delete $atom_y{$i};
	delete $atom_z{$i};
	delete $atom_occ{$i};
	delete $footnote_number{$i};
	delete $B_or_U{$i};
	if ($sigatm_flag > 0 && ($atom_or_het{$i} ne 'TER' ||

	  ($print_ter ne 'no' && $print_ter ne 'comment'))) {
	    $X = $sigatm_point{($atom_number{$i} . '|' . $cur_model)};
	    $ys = '';
	    if (substr($sig_y{$X}, 1, 1) ne ' ') {
		$ys = "\n                                     ";
	    }
	    $zs = '';
	    if (substr($sig_z{$X}, 1, 1) ne ' ') {
		$zs = "\n                                             ";
	    }
	    $os = '';
	    if (substr($sig_occ{$X}, 1, 1) ne ' ') {
		$os =

		  "\n                                                     ";
	    }
	    $bs = '';
	    if (substr($sig_temp{$X}, 1, 1) ne ' ') {
		$bs =

		  "\n                                                           ";
	    }
	    if ($X ne '') {	#???
		printf

		  "                             %8s%s%8s%s%8s%s%6s%s%6s\n",

		  $sig_x{$X}, $ys, $sig_y{$X}, $zs, $sig_z{$X}, $os,

		  $sig_occ{$X}, $bs, $sig_temp{$X};
	    }
	    else {
		printf

		  (("                                 .        .        .       .      .\n"));
	    }
	}
	# end if (sigatm_flag > 0)
	if ($aniso_flag > 0 && ($atom_or_het{$i} ne 'TER' ||

	  ($print_ter ne 'no' && $print_ter ne 'comment'))) {
	    $X = $aniso_point{($atom_number{$i} . '|' . $cur_model)};
	    if ($X ne '') {	#???
		printf "                       %7g %7g %7g %7g %7g %7g \n",

		  $atom_U11{$X} / 10000, $atom_U12{$X} / 10000,

		  $atom_U13{$X} / 10000, $atom_U22{$X} / 10000,

		  $atom_U23{$X} / 10000, $atom_U33{$X} / 10000;
	    }
	    else {
		printf (("    .       .       .       .       .       .\n"));
	    }
	}
	# end if (aniso_flag > 0)
	if ($siguij_flag > 0 && ($atom_or_het{$i} ne 'TER' ||

	  ($print_ter ne 'no' && $print_ter ne 'comment'))) {
	    $X = $siguij_point{($atom_number{$i} . '|' . $cur_model)};
	    if ($X ne '') {	#???
		printf "                       %7g %7g %7g %7g %7g %7g \n",

		  $sig_U11{$X} / 10000, $sig_U12{$X} / 10000,

		  $sig_U13{$X} / 10000, $sig_U22{$X} / 10000,

		  $sig_U23{$X} / 10000, $sig_U33{$X} / 10000;
	    }
	    else {
		printf (("    .       .       .       .       .       .\n"));
	    }
	}
	# end if (siguij_flag > 0)

	;
    }
}
# end for (i=1; i < atom_flag; ++i)
if ($nmr_flag > 0 || $vol_flag > 0) {
    if ($b_warn > 0 || $occ_warn > 0) {
	$warning_list{++$warning_flag} =

	  "#=# ATOM_SITE: Unusual occupancies or B-values, read REMARKs\n";
    }
}
else {
    if ($occ_warn > 0) {
	$warning_list{++$warning_flag} =

	  "#=# ATOM_SITE: Unusual occupancies, read REMARKs\n";
    }
}
#
#
#
# Process HEADER records for DATABASE_2
{
    printf (("\n\n\n"));
    printf (("####################\n"));
    printf (("#                  #\n"));
    printf (("# DATABASE_2       #\n"));
    printf (("#                  #\n"));
    printf (("####################\n\n\n"));
    printf (("_database_2.database_id			PDB \n"));
    printf "_database_2.database_code		%4s \n", $head_PDB_code;
}
#
#
#
# Process REVDAT records results
{
    if ($revdat_flag ne '1' || $s_o_flag != 0) {	#???
	printf (("\n\n\n"));
	printf (("####################\n"));
	printf (("#                  #\n"));
	printf (("# DATABASE_PDB_REV #\n"));
	printf (("#                  #\n"));
	printf (("####################\n\n\n"));
	{
	    $null = (@odmy = split(/-/, $head_dep_date, 9999));
	    if ($null >= 1 && $odmy[$null] eq '' && '-' eq ' ') {
		--$null;
	    }
	}
	printf (("\nloop_\n"));
	printf (("_database_PDB_rev.date_original\n"));
	printf (("_database_PDB_rev.num\n"));
	printf (("_database_PDB_rev.date\n"));
	printf (("_database_PDB_rev.mod_type\n"));
	printf (("_database_PDB_rev.status\n"));
	printf (("_database_PDB_rev.replaced_by\n"));
	printf (("_database_PDB_rev.replaces\n"));

	# Process SPRSDE and OBSLTE records

	for ($i = 1; $i <= $s_o_flag; ++$i) {
	    $my_date = $s_o_date{$i};
	    $my_num = $rev_mod_of_date{$my_date};
	    if ($my_num eq '' || $my_num eq ' ') {	#???	#???
		++$revdat_flag;
		$rev_mod_number{$revdat_flag - 1} = '?';
		$my_num = '?';
		if ($my_date eq '         ') {
		    $my_date = '??-???-??';
		}
		{
		    $cal_date = (@dmy = split(/-/, $my_date, 9999));
		    if ($cal_date >= 1 && $dmy[$cal_date] eq '' &&

		      '-' eq ' ') {
			--$cal_date;
		    }
		}
		$rev_date_year{$revdat_flag - 1} = $yyyy{$dmy[3] + 0};
		$rev_date_mon{$revdat_flag - 1} = $mmm2mm{$dmy[2]};
		$rev_date_day{$revdat_flag - 1} = $dmy[1];
		$rev_cont_flag{$revdat_flag - 1} = '  ';
		$rev_name{$revdat_flag - 1} = '?';
		$rev_type{$revdat_flag - 1} = '?';
		$rev_rec_corr{$revdat_flag - 1} = ' ';
		$warning_list{++$warning_flag} =

		  "#=# DATABASE_PDB_REV: Unidentified date from OBSLTE or SPRSDE\n";
	    }
	    $rev_s_o{$my_num} = ($rev_s_o{$my_num} . ' ' . $i);
	}
	for ($j = 1; $j < $revdat_flag; ++$j) {
	    if ($rev_cont_flag{$j} eq '  ') {
		$my_num = $rev_mod_number{$j};
		if ($my_num ne '?') {	#???
		    $xrev_s_o = $rev_s_o{$my_num + 0};
		}
		else {
		    $xrev_s_o = $rev_s_o{$my_num};
		}
		{
		    $my_s_o = (@x_s_o = split(' ', $xrev_s_o, 9999));
		    if ($my_s_o >= 1 && $x_s_o[$my_s_o] eq '' && ' ' eq ' ') {
			--$my_s_o;
		    }
		}
		if ($my_s_o == 0 || $x_s_o[1] eq '' || $x_s_o[1] eq ' ') {
		    printf ' %4s-%2s-%2s ', $yyyy{$odmy[3] + 0},

		      $mmm2mm{$odmy[2]}, $odmy[1];
		    printf "%3s  %4s-%2s-%2s  %1s . . .\n",

		      $rev_mod_number{$j}, $rev_date_year{$j},

		      $rev_date_mon{$j}, $rev_date_day{$j}, $rev_type{$j};
		}
		else {
		    for ($k = 1; $k <= $my_s_o; ++$k) {
			$i = $x_s_o[$k];
			$s_o_s = '.';
			$s_o_o = '.';
			$s_o_stat = '.';
			$my_type = $s_o_type{$i};
			$my_list = $s_o_list{$i};
			{
			    $my_count = (@my_revs = split(' ', $my_list,

			      9999));
			    if ($my_count >= 1 && $my_revs[$my_count] eq '' &&

			      ' ' eq ' ') {
				--$my_count;
			    }
			}
			for ($l = 1; $l <= $my_count; ++$l) {
			    if ($my_type eq 'OBSLTE') {
				$s_o_o = $my_revs[$l];
				$s_o_stat = 'obsolete';
			    }
			    else {
				$s_o_s = $my_revs[$l];
			    }
			    printf ' %4s-%2s-%2s ', $yyyy{$odmy[3] + 0},

			      $mmm2mm{$odmy[2]}, $odmy[1];
			    printf "%3s  %4s-%2s-%2s  %1s %s %s %s\n",

			      $rev_mod_number{$j}, $rev_date_year{$j},

			      $rev_date_mon{$j}, $rev_date_day{$j},

			      $rev_type{$j}, $s_o_stat, $s_o_o, $s_o_s;
			}
		    }
		    # end for (k = 1; k <= my_s_o; ++k)

		    ;
		}
		# end (my_s_o == 0 || x_s_o[1] == "" || x_s_o[1] == " ")

		;
	    }
	    # end if ( rev_cont_flag[j] == " " )

	    ;
	}
	$pdb_reloop = 0;
	for ($j = 1; $j < $revdat_flag; ++$j) {
	    # Additional data items if more than initial entry
	    {
		$rev_num_change = (@rev_recs = split(' ', $rev_rec_corr{$j},

		  9999));
		if ($rev_num_change >= 1 &&

		  $rev_recs[$rev_num_change] eq '' && ' ' eq ' ') {
		    --$rev_num_change;
		}
	    }
	    for ($k = 1; $k <= $rev_num_change; ++$k) {
		if ($pdb_reloop == 0) {
		    printf (("\nloop_\n"));
		    printf (("_database_PDB_rev_record.rev_num\n"));
		    printf (("_database_PDB_rev_record.details\n"));
		    printf (("_database_PDB_rev_record.type\n"));
		    $pdb_reloop++;
		}
		# end if (pdb_reloop == 0 )
		printf " %3s        %5s        %-8s\n", $rev_mod_number{$j},

		  $rev_name{$j}, $rev_recs[$k];
	    }
	}
	# end for (j=1; j < revdat_flag; ++j) 

	;
    }
    # end if (revdat_flag != "1")
    #
    #  Display geometry data items here

    if ($verbose eq 'yes') {
	printf (("\n_geom.details           ? \n"));

	printf (("\nloop_\n"));
	printf (("_geom_angle.atom_site_id_1 \n"));
	printf (("_geom_angle.atom_site_id_2 \n"));
	printf (("_geom_angle.atom_site_id_3 \n"));
	printf (("_geom_angle.value \n"));
	printf (("_geom_angle.site_symmetry_1 \n"));
	printf (("_geom_angle.site_symmetry_2 \n"));
	printf (("_geom_angle.site_symmetry_3 \n"));
	printf (("_geom_angle.publ_flag \n"));
	if ($model_flags eq 'yes') {
	    printf ((('_geom_angle.' . $pdb2cif_prefix .

	      "auth_id_in_model_1 \n")));
	    printf ((('_geom_angle.' . $pdb2cif_prefix .

	      "auth_id_in_model_2 \n")));
	    printf ((('_geom_angle.' . $pdb2cif_prefix .

	      "auth_id_in_model_3 \n")));
	    printf ((('_geom_angle.' . $pdb2cif_prefix .

	      "auth_model_id \n")));
	    printf ((('_geom_angle.' . $pdb2cif_prefix .

	      "id_in_model_1 \n")));
	    printf ((('_geom_angle.' . $pdb2cif_prefix .

	      "id_in_model_2 \n")));
	    printf ((('_geom_angle.' . $pdb2cif_prefix .

	      "id_in_model_3 \n")));
	    printf ((('_geom_angle.' . $pdb2cif_prefix .

	      "label_model_id \n")));
	}
	printf ((" ?  ?  ?  ?  ?  ?  ?  ? \n"));
	if ($model_flags eq 'yes') {
	    printf (("   ?  ?  ?  ?  ?  ?  ?  ? \n"));
	}

	printf (("\nloop_\n"));
	printf (("_geom_bond.atom_site_id_1 \n"));
	printf (("_geom_bond.atom_site_id_2 \n"));
	printf (("_geom_bond.dist \n"));
	printf (("_geom_bond.site_symmetry_1 \n"));
	printf (("_geom_bond.site_symmetry_2 \n"));
	printf (("_geom_bond.publ_flag \n"));
	if ($model_flags eq 'yes') {
	    printf ((('_geom_bond.' . $pdb2cif_prefix .

	      "auth_id_in_model_1 \n")));
	    printf ((('_geom_bond.' . $pdb2cif_prefix .

	      "auth_id_in_model_2 \n")));
	    printf ((('_geom_bond.' . $pdb2cif_prefix . "auth_model_id \n")));
	    printf ((('_geom_bond.' . $pdb2cif_prefix . "id_in_model_1 \n")));
	    printf ((('_geom_bond.' . $pdb2cif_prefix . "id_in_model_2 \n")));
	    printf ((('_geom_bond.' . $pdb2cif_prefix .

	      "label_model_id \n")));
	}
	printf ((" ?  ?  ?  ?  ?  ?  \n"));
	if ($model_flags eq 'yes') {
	    printf (("   ?  ?  ?  ?  ?  ? \n"));
	}

	printf (("\nloop_\n"));
	printf (("_geom_contact.atom_site_id_1 \n"));
	printf (("_geom_contact.atom_site_id_2 \n"));
	printf (("_geom_contact.dist \n"));
	printf (("_geom_contact.publ_flag \n"));
	printf (("_geom_contact.site_symmetry_1 \n"));
	printf (("_geom_contact.site_symmetry_2 \n"));
	if ($model_flags eq 'yes') {
	    printf ((('_geom_contact.' . $pdb2cif_prefix .

	      "auth_id_in_model_1 \n")));
	    printf ((('_geom_contact.' . $pdb2cif_prefix .

	      "auth_id_in_model_2 \n")));
	    printf ((('_geom_contact.' . $pdb2cif_prefix .

	      "auth_model_id \n")));
	    printf ((('_geom_contact.' . $pdb2cif_prefix .

	      "id_in_model_1 \n")));
	    printf ((('_geom_contact.' . $pdb2cif_prefix .

	      "id_in_model_2 \n")));
	    printf ((('_geom_contact.' . $pdb2cif_prefix .

	      "label_model_id \n")));
	}
	printf ((" ?  ?  ?  ?  ?  ?  \n"));
	if ($model_flags eq 'yes') {
	    printf (("   ?  ?  ?  ?  ?  ? \n"));
	}

	printf (("\nloop_\n"));
	printf (("_geom_torsion.atom_site_id_1 \n"));
	printf (("_geom_torsion.atom_site_id_2 \n"));
	printf (("_geom_torsion.atom_site_id_3 \n"));
	printf (("_geom_torsion.atom_site_id_4 \n"));
	printf (("_geom_torsion.value \n"));
	printf (("_geom_torsion.publ_flag \n"));
	printf (("_geom_torsion.site_symmetry_1 \n"));
	printf (("_geom_torsion.site_symmetry_2 \n"));
	printf (("_geom_torsion.site_symmetry_3 \n"));
	printf (("_geom_torsion.site_symmetry_4 \n"));
	if ($model_flags eq 'yes') {
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "auth_id_in_model_1 \n")));
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "auth_id_in_model_2 \n")));
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "auth_id_in_model_3 \n")));
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "auth_id_in_model_4 \n")));
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "auth_model_id \n")));
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "id_in_model_1 \n")));
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "id_in_model_2 \n")));
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "id_in_model_3 \n")));
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "id_in_model_4 \n")));
	    printf ((('_geom_torsion.' . $pdb2cif_prefix .

	      "label_model_id \n")));
	}

	printf ((" ?  ?  ?  ?  ?  ?  ?  ?  ?  ? \n"));
	if ($model_flags eq 'yes') {
	    printf (("   ?  ?  ?  ?  ?  ?  ?  ?  ?  ? \n"));
	}
    }

    #
    #
    #
    # Process HEADER records for STRUCT_BIOL
    {
	printf (("\n\n\n"));
	printf (("####################\n"));
	printf (("#                  #\n"));
	printf (("# STRUCT_BIOL      #\n"));
	printf (("#                  #\n"));
	printf (("####################\n\n\n"));

	printf (("\nloop_\n"));
	printf (("_struct_biol.id\n"));
	printf (("_struct_biol.details\n"));
	printf " %4s '%s' \n", $head_PDB_code, $head_funct_class;
	for ($i = 1; $i <= $num_poly_ents; ++$i) {
	    printf " %4s   . \n", $entities{$i};
	}

	printf (("\nloop_\n"));
	printf (("_struct_biol_gen.biol_id\n"));
	printf (("_struct_biol_gen.asym_id\n"));
	printf (("_struct_biol_gen.symmetry\n"));
	printf (("_struct_biol_gen.details\n"));
	for ($i = 1; $i < $num_ents; ++$i) {
	    printf " %4s   %4s   1_555   .\n", $head_PDB_code, $entities{$i};
	}
	for ($i = 1; $i <= $num_poly_ents; ++$i) {
	    printf " %4s   %4s   1_555   .\n", $entities{$i}, $entities{$i};
	}
    }
    #
    #  process pending connect information
    #
    {
	if ($connect_flag + 0 > 0 || $ssbond_flag > 1) {
	    printf (("\n\n\n"));
	    printf (("##############################\n"));
	    printf (("#                            #\n"));
	    printf (("# STRUCT_CONN_TYPE           #\n"));
	    printf (("#                            #\n"));
	    printf (("##############################\n\n\n"));

	    printf (("\nloop_\n"));
	    printf (("_struct_conn_type.id\n"));
	    printf (("_struct_conn_type.criteria\n"));
	    printf (("_struct_conn_type.reference\n"));
	    printf (("  .     'unknown bond type from PDB entry' ?\n"));
	    printf (("saltbr  'salt bridge from PDB entry'       ?\n"));
	    printf (("hydrog  'hydrogen bond from PDB entry'     ?\n"));
	    if ($ssbond_flag > 1) {
		printf (("disulf  'disulfide bridge from PDB entry'  ?\n"));
	    }

	    printf (("\n\n\n"));
	    printf (("##############################\n"));
	    printf (("#                            #\n"));
	    printf (("# STRUCT_CONN                #\n"));
	    printf (("#                            #\n"));
	    printf (("##############################\n\n\n"));

	    printf (("\nloop_\n"));
	    printf (("_struct_conn.id\n"));
	    printf (("_struct_conn.conn_type_id\n"));
	    printf (("_struct_conn.ptnr1_label_comp_id\n"));
	    printf (("_struct_conn.ptnr1_label_asym_id\n"));
	    printf (("_struct_conn.ptnr1_auth_seq_id\n"));
	    printf (("_struct_conn.ptnr1_label_atom_id\n"));
	    printf (("_struct_conn.ptnr1_label_alt_id\n"));
	    printf (("_struct_conn.ptnr1_role\n"));
	    printf (("_struct_conn.ptnr1_symmetry\n"));
	    printf (("_struct_conn.ptnr2_label_comp_id\n"));
	    printf (("_struct_conn.ptnr2_label_asym_id\n"));
	    printf (("_struct_conn.ptnr2_auth_seq_id\n"));
	    printf (("_struct_conn.ptnr2_label_atom_id\n"));
	    printf (("_struct_conn.ptnr2_label_alt_id\n"));
	    printf (("_struct_conn.ptnr2_role\n"));
	    printf (("_struct_conn.ptnr2_symmetry\n"));
	    printf ((('_struct_conn.' . $pdb2cif_prefix .

	      "ptnr1_atom_site_id\n")));
	    printf (("_struct_conn.ptnr1_label_seq_id\n"));
	    printf ((('_struct_conn.' . $pdb2cif_prefix .

	      "ptnr2_atom_site_id\n")));
	    printf (("_struct_conn.ptnr2_label_seq_id\n"));
	    if ($model_flags eq 'yes') {
		printf ((('_struct_conn.' . $pdb2cif_prefix .

		  "ptnr1_id_in_model\n")));
		printf ((('_struct_conn.' . $pdb2cif_prefix .

		  "ptnr1_auth_id_in_model\n")));
		printf ((('_struct_conn.' . $pdb2cif_prefix .

		  "ptnr2_id_in_model\n")));
		printf ((('_struct_conn.' . $pdb2cif_prefix .

		  "ptnr2_auth_id_in_model\n")));
		printf ((('_struct_conn.' . $pdb2cif_prefix .

		  "label_model_id\n")));
	    }
	}

	# parse field

	$xconnect_flag = 0;
	for ($icon = 1; $icon <= $connect_flag; ++$icon) {
	    $connect_source = substr($connect_save{$icon}, 7, 5);
	    $itarget = 0;
	    for ($ipos = 12; $ipos <= 57; $ipos += 5) {
		$connect_tar{++$itarget} = substr($connect_save{$icon}, $ipos,

		  5);
		$c1 = 0 + $connect_source;
		$c2 = 0 + $connect_tar{$itarget};
		if ($connect_tar{$itarget} ne '     ') {
		    if (($chain_id{$atom_point{$c1}} .

		      $residue_seq_number{$atom_point{$c1}} .

		      $atom_name{$atom_point{$c1}} .

		      $atom_alt_location{$atom_point{$c1}}) gt	#???

		      ($chain_id{$atom_point{$c2}} .

		      $residue_seq_number{$atom_point{$c2}} .

		      $atom_name{$atom_point{$c2}} .

		      $atom_alt_location{$atom_point{$c2}})) {
			$ct = $c1;
			$c1 = $c2;
			$c2 = $ct;
		    }
		    # end if ((chain_id[atom_point[c1]] residue_seq_number[atom_point[c1]] ...
		    $conect_ptnr1_label_res_id =

		      $residue_name{$atom_point{$c1}};
		    $conect_ptnr1_label_asym_id = $chain_id{$atom_point{$c1}};
		    $conect_ptnr1_auth_seq_id =

		      $residue_seq_number{$atom_point{$c1}};
		    $conect_ptnr1_label_atom_id =

		      $atom_name{$atom_point{$c1}};
		    $conect_ptnr1_label_alt_id =

		      $atom_alt_location{$atom_point{$c1}};
		    $conect_ptnr2_label_res_id =

		      $residue_name{$atom_point{$c2}};
		    $conect_ptnr2_label_asym_id = $chain_id{$atom_point{$c2}};
		    $conect_ptnr2_auth_seq_id =

		      $residue_seq_number{$atom_point{$c2}};
		    $conect_ptnr2_label_atom_id =

		      $atom_name{$atom_point{$c2}};
		    $conect_ptnr2_label_alt_id =

		      $atom_alt_location{$atom_point{$c2}};
		    $connect_bond = ($conect_ptnr1_label_res_id .

		      $conect_ptnr1_label_asym_id . $conect_ptnr1_auth_seq_id

		      . $conect_ptnr1_label_atom_id .

		      $conect_ptnr1_label_alt_id . ' ' .

		      $conect_ptnr2_label_res_id . $conect_ptnr2_label_asym_id

		      . $conect_ptnr2_auth_seq_id .

		      $conect_ptnr2_label_atom_id .

		      $conect_ptnr2_label_alt_id);
		    if ($bond_count{$connect_bond} + 0 == 0) {
			$aslist{$connect_bond} = '|';
			++$bond_count{$connect_bond};
			$conect_xtype{$connect_bond} =

			  $connect_types[$itarget];
			$conect_symm_1{$connect_bond} = '   .   ';
			$conect_symm_2{$connect_bond} = '   .   ';
			$conect_role_1{$connect_bond} = '.';
			$conect_role_2{$connect_bond} = '.';
			$conect_epsn1{$connect_bond} =

			  $ent_pol_seq_num{($conect_ptnr1_label_asym_id .

			  $conect_ptnr1_auth_seq_id)};
			$conect_epsn2{$connect_bond} =

			  $ent_pol_seq_num{($conect_ptnr2_label_asym_id .

			  $conect_ptnr2_auth_seq_id)};
			if ($conect_epsn1{$connect_bond} eq '') {
			    $conect_epsn1{$connect_bond} = '. ';
			}
			if ($conect_epsn2{$connect_bond} eq '') {
			    $conect_epsn2{$connect_bond} = '. ';
			}
			$conect_c1_label_res_id{$connect_bond} =

			  $conect_ptnr1_label_res_id;
			$conect_c1_label_asym_id{$connect_bond} =

			  $conect_ptnr1_label_asym_id;
			$conect_c1_auth_seq_id{$connect_bond} =

			  $conect_ptnr1_auth_seq_id;
			$conect_c1_label_atom_id{$connect_bond} =

			  $conect_ptnr1_label_atom_id;
			$conect_c1_label_alt_id{$connect_bond} =

			  $conect_ptnr1_label_alt_id;
			$conect_c2_label_res_id{$connect_bond} =

			  $conect_ptnr2_label_res_id;
			$conect_c2_label_asym_id{$connect_bond} =

			  $conect_ptnr2_label_asym_id;
			$conect_c2_auth_seq_id{$connect_bond} =

			  $conect_ptnr2_auth_seq_id;
			$conect_c2_label_atom_id{$connect_bond} =

			  $conect_ptnr2_label_atom_id;
			$conect_c2_label_alt_id{$connect_bond} =

			  $conect_ptnr2_label_alt_id;
			$next_bond{++$xconnect_flag} = $connect_bond;
		    }
		    # end if (bond_count[connect_bond]+0 == 0)
		    if (index($aslist{$connect_bond},

		      ('|' . $c1 . ' ' . $c2 . '|')) == 0) {
			$aslist{$connect_bond} = ($aslist{$connect_bond} . $c1

			  . ' ' . $c2 . '|');
		    }
		}
		# end if (c2 != "     ")

		;
	    }
	    # end for (ipos = 12 ; ipos <= 57; ipos += 5)

	    ;
	}
	# end for(icon = 1; icon <= connect_flag; ++ icon)
	#
	#  add disulphides if present to CONECT list
	#

	if ($ssbond_flag > 1) {
	    for ($i = 1; $i < $ssbond_flag; ++$i) {
		$ssbs1 = $ssbond_symm_1{$i};
		$ssbs2 = $ssbond_symm_2{$i};
		$connect_bond = ($ssbond_res_name_beg{$i} .

		  $ssbond_chain_id_beg{$i} . $ssbond_res_seq_num_beg{$i} .

		  ' SG ' . '.' . ' ' . $ssbond_res_name_end{$i} .

		  $ssbond_chain_id_end{$i} . $ssbond_res_seq_num_end{$i} .

		  ' SG ' . '.');
		if (($ssbond_chain_id_beg{$i} . $ssbond_res_seq_num_beg{$i} .

		  ' SG ' . '.') gt ($ssbond_chain_id_end{$i} .

		  $ssbond_res_seq_num_end{$i} . ' SG ' . '.')) {
		    $ssbs1 = $ssbond_symm_2{$i};
		    $ssbs2 = $ssbond_symm_1{$i};
		    $connect_bond = ($ssbond_res_name_end{$i} .

		      $ssbond_chain_id_end{$i} . $ssbond_res_seq_num_end{$i} .

		      ' SG ' . '.' . ' ' . $ssbond_res_name_beg{$i} .

		      $ssbond_chain_id_beg{$i} . $ssbond_res_seq_num_beg{$i} .

		      ' SG ' . '.');
		}
		# end if ((ssbond_chain_id_beg[i] ssbond_res_seq_num_beg[i] ...
		if ($bond_count{$connect_bond} + 0 == 0) {
		    $warning_list{++$warning_flag} =

		      ('#=# STRUCT_CONN: No connect for disulf ' .

		      $connect_bond . " \n");
		    $aslist{$connect_bond} = ('|. .|');
		    ++$bond_count{$connect_bond};
		    $conect_xtype{$connect_bond} = 'disulf';
		    $conect_symm_1{$connect_bond} = $ssbs1;
		    $conect_symm_2{$connect_bond} = $ssbs2;
		    $conect_ptnr1_label_asym_id = substr($connect_bond, 4, 1);
		    $conect_ptnr1_auth_seq_id = substr($connect_bond, 5, 5);
		    $conect_ptnr2_label_asym_id = substr($connect_bond, 19,

		      1);
		    $conect_ptnr2_auth_seq_id = substr($connect_bond, 20, 5);
		    $conect_role_1{$connect_bond} = '.';
		    $conect_role_2{$connect_bond} = '.';
		    $conect_epsn1{$connect_bond} =

		      $ent_pol_seq_num{($conect_ptnr1_label_asym_id .

		      $conect_ptnr1_auth_seq_id)};
		    $conect_epsn2{$connect_bond} =

		      $ent_pol_seq_num{($conect_ptnr2_label_asym_id .

		      $conect_ptnr2_auth_seq_id)};
		    if ($conect_epsn1{$connect_bond} eq '') {
			$conect_epsn1{$connect_bond} = '. ';
		    }
		    if ($conect_epsn2{$connect_bond} eq '') {
			$conect_epsn2{$connect_bond} = '. ';
		    }
		    $conect_c1_label_res_id{$connect_bond} =

		      substr($connect_bond, 1, 3);
		    $conect_c1_label_asym_id{$connect_bond} =

		      $conect_ptnr1_label_asym_id;
		    $conect_c1_auth_seq_id{$connect_bond} =

		      $conect_ptnr1_auth_seq_id;
		    $conect_c1_label_atom_id{$connect_bond} =

		      substr($connect_bond, 10, 4);
		    $conect_c1_label_alt_id{$connect_bond} =

		      substr($connect_bond, 14, 1);
		    $conect_c2_label_res_id{$connect_bond} =

		      substr($connect_bond, 16, 3);
		    $conect_c2_label_asym_id{$connect_bond} =

		      $conect_ptnr2_label_asym_id;
		    $conect_c2_auth_seq_id{$connect_bond} =

		      $conect_ptnr2_auth_seq_id;
		    $conect_c2_label_atom_id{$connect_bond} =

		      substr($connect_bond, 25, 4);
		    $conect_c2_label_alt_id{$connect_bond} =

		      substr($connect_bond, 29, 1);
		    $next_bond{++$xconnect_flag} = $connect_bond;
		}
		else {
		    $conect_xtype{$connect_bond} = 'disulf';
		    $conect_symm_1{$connect_bond} = $ssbs1;
		    $conect_symm_2{$connect_bond} = $ssbs2;
		}
		# end if (bond_count[connect_bond]+0 == 0)

		;
	    }
	    # end for (i=1; i < ssbond_flag; ++i)

	    ;
	}
	# end if (ssbond_flag > 1 )
	if ($sltbrg_flag > 1) {
	    for ($i = 1; $i < $sltbrg_flag; ++$i) {
		$sbra1 = $sb_atom_beg{$i};
		$sbra2 = $sb_atom_end{$i};
		$sbrs1 = $sb_symm_1{$i};
		$sbrs2 = $sb_symm_2{$i};
		$sbral1 = $sb_alt_loc_beg{$i};
		$sbral2 = $sb_alt_loc_end{$i};
		$connect_bond = ($sb_res_name_beg{$i} . $sb_chain_id_beg{$i} .

		  $sb_res_seq_num_beg{$i} . $sbra1 . $sbral1 . ' ' .

		  $sb_res_name_end{$i} . $sb_chain_id_end{$i} .

		  $sb_res_seq_num_end{$i} . $sbra2 . $sbral2);
		if (($sb_chain_id_beg{$i} . $sb_res_seq_num_beg{$i} . $sbra1 .

		  $sbral1) gt ($sb_chain_id_end{$i} . $sb_res_seq_num_end{$i}	#???

		  . $sbra2 . $sbral2)) {
		    $sbra1 = $sb_atom_end{$i};
		    $sbra2 = $sb_atom_beg{$i};
		    $sbrs1 = $sb_symm_2{$i};
		    $sbrs2 = $sb_symm_1{$i};
		    $sbral1 = $sb_alt_loc_end{$i};
		    $sbral2 = $sb_alt_loc_beg{$i};
		    $connect_bond = ($sb_res_name_end{$i} .

		      $sb_chain_id_end{$i} . $sb_res_seq_num_end{$i} . $sbra1

		      . $sbral1 . ' ' . $sb_res_name_beg{$i} .

		      $sb_chain_id_beg{$i} . $sb_res_seq_num_beg{$i} . $sbra2

		      . $sbral2);
		}
		# end if ((sb_chain_id_beg[i] sb_res_seq_num_beg[i] ...
		if ($bond_count{$connect_bond} + 0 == 0) {
		    if ($sbrs1 eq $sbrs2) {	#???
			$warning_list{++$warning_flag} =

			  ('#=# STRUCT_CONN: No connect for saltbr ' .

			  $connect_bond . " \n");
		    }
		    $aslist{$connect_bond} = ('|. .|');
		    ++$bond_count{$connect_bond};
		    $conect_xtype{$connect_bond} = 'saltbr';
		    $conect_symm_1{$connect_bond} = $sbrs1;
		    $conect_symm_2{$connect_bond} = $sbrs2;
		    $conect_ptnr1_label_asym_id = substr($connect_bond, 4, 1);
		    $conect_ptnr1_auth_seq_id = substr($connect_bond, 5, 5);
		    $conect_ptnr2_label_asym_id = substr($connect_bond, 19,

		      1);
		    $conect_ptnr2_auth_seq_id = substr($connect_bond, 20, 5);
		    $conect_role_1{$connect_bond} = '.';
		    $conect_role_2{$connect_bond} = '.';
		    $conect_epsn1{$connect_bond} =

		      $ent_pol_seq_num{($conect_ptnr1_label_asym_id .

		      $conect_ptnr1_auth_seq_id)};
		    $conect_epsn2{$connect_bond} =

		      $ent_pol_seq_num{($conect_ptnr2_label_asym_id .

		      $conect_ptnr2_auth_seq_id)};
		    if ($conect_epsn1{$connect_bond} eq '') {
			$conect_epsn1{$connect_bond} = '. ';
		    }
		    if ($conect_epsn2{$connect_bond} eq '') {
			$conect_epsn2{$connect_bond} = '. ';
		    }
		    $conect_c1_label_res_id{$connect_bond} =

		      substr($connect_bond, 1, 3);
		    $conect_c1_label_asym_id{$connect_bond} =

		      $conect_ptnr1_label_asym_id;
		    $conect_c1_auth_seq_id{$connect_bond} =

		      $conect_ptnr1_auth_seq_id;
		    $conect_c1_label_atom_id{$connect_bond} =

		      substr($connect_bond, 10, 4);
		    $conect_c1_label_alt_id{$connect_bond} =

		      substr($connect_bond, 14, 1);
		    $conect_c2_label_res_id{$connect_bond} =

		      substr($connect_bond, 16, 3);
		    $conect_c2_label_asym_id{$connect_bond} =

		      $conect_ptnr2_label_asym_id;
		    $conect_c2_auth_seq_id{$connect_bond} =

		      $conect_ptnr2_auth_seq_id;
		    $conect_c2_label_atom_id{$connect_bond} =

		      substr($connect_bond, 25, 4);
		    $conect_c2_label_alt_id{$connect_bond} =

		      substr($connect_bond, 29, 1);
		    $next_bond{++$xconnect_flag} = $connect_bond;
		}
		else {
		    if ($conect_xtype{$connect_bond} ne 'saltbr') {
			$warning_list{++$warning_flag} =

			  ('#=# STRUCT_CONN: Changed bond ' . $connect_bond .

			  " to saltbr\n");
		    }
		    $conect_xtype{$connect_bond} = 'saltbr';
		    $conect_symm_1{$connect_bond} = $sbrs1;
		    $conect_symm_2{$connect_bond} = $sbrs2;
		}
		# end if (bond_count[connect_bond]+0 == 0)

		;
	    }
	    # end for (i=1; i < sltbrg_flag; ++i)

	    ;
	}
	# end if (sltbrg_flag > 1 )
	if ($link_flag > 1) {
	    for ($i = 1; $i < $link_flag; ++$i) {
		$lnka1 = $lk_atom_beg{$i};
		$lnka2 = $lk_atom_end{$i};
		$lnks1 = $lk_symm_1{$i};
		$lnks2 = $lk_symm_2{$i};
		$lnkal1 = $lk_alt_loc_beg{$i};
		$lnkal2 = $lk_alt_loc_end{$i};
		$connect_bond = ($lk_res_name_beg{$i} . $lk_chain_id_beg{$i} .

		  $lk_res_seq_num_beg{$i} . $lnka1 . $lnkal1 . ' ' .

		  $lk_res_name_end{$i} . $lk_chain_id_end{$i} .

		  $lk_res_seq_num_end{$i} . $lnka2 . $lnkal2);
		if (($lk_chain_id_beg{$i} . $lk_res_seq_num_beg{$i} . $lnka1 .

		  $lnkal1) gt ($lk_chain_id_end{$i} . $lk_res_seq_num_end{$i}	#???

		  . $lnka2 . $lnkal2)) {
		    $lnka1 = $lk_atom_end{$i};
		    $lnka2 = $lk_atom_beg{$i};
		    $lnks1 = $lk_symm_2{$i};
		    $lnks2 = $lk_symm_1{$i};
		    $lnkal1 = $lk_alt_loc_end{$i};
		    $lnkal2 = $lk_alt_loc_beg{$i};
		    $connect_bond = ($lk_res_name_end{$i} .

		      $lk_chain_id_end{$i} . $lk_res_seq_num_end{$i} . $lnka1

		      . $lnkal1 . ' ' . $lk_res_name_beg{$i} .

		      $lk_chain_id_beg{$i} . $lk_res_seq_num_beg{$i} . $lnka2

		      . $lnkal2);
		}
		# end if ((lk_chain_id_beg[i] lk_res_seq_num_beg[i] ...
		if ($bond_count{$connect_bond} + 0 == 0) {
		    if ($lnks1 eq $lnks2) {	#???
			$warning_list{++$warning_flag} =

			  ('#=# STRUCT_CONN: No connect for link ' .

			  $connect_bond . " \n");
		    }
		    $aslist{$connect_bond} = ('|. .|');
		    ++$bond_count{$connect_bond};
		    $conect_xtype{$connect_bond} = '.';
		    $conect_symm_1{$connect_bond} = $lnks1;
		    $conect_symm_2{$connect_bond} = $lnks2;
		    $conect_ptnr1_label_asym_id = substr($connect_bond, 4, 1);
		    $conect_ptnr1_auth_seq_id = substr($connect_bond, 5, 5);
		    $conect_ptnr2_label_asym_id = substr($connect_bond, 19,

		      1);
		    $conect_ptnr2_auth_seq_id = substr($connect_bond, 20, 5);
		    $conect_role_1{$connect_bond} = '.';
		    $conect_role_2{$connect_bond} = '.';
		    $conect_epsn1{$connect_bond} =

		      $ent_pol_seq_num{($conect_ptnr1_label_asym_id .

		      $conect_ptnr1_auth_seq_id)};
		    $conect_epsn2{$connect_bond} =

		      $ent_pol_seq_num{($conect_ptnr2_label_asym_id .

		      $conect_ptnr2_auth_seq_id)};
		    if ($conect_epsn1{$connect_bond} eq '') {
			$conect_epsn1{$connect_bond} = '. ';
		    }
		    if ($conect_epsn2{$connect_bond} eq '') {
			$conect_epsn2{$connect_bond} = '. ';
		    }
		    $conect_c1_label_res_id{$connect_bond} =

		      substr($connect_bond, 1, 3);
		    $conect_c1_label_asym_id{$connect_bond} =

		      $conect_ptnr1_label_asym_id;
		    $conect_c1_auth_seq_id{$connect_bond} =

		      $conect_ptnr1_auth_seq_id;
		    $conect_c1_label_atom_id{$connect_bond} =

		      substr($connect_bond, 10, 4);
		    $conect_c1_label_alt_id{$connect_bond} =

		      substr($connect_bond, 14, 1);
		    $conect_c2_label_res_id{$connect_bond} =

		      substr($connect_bond, 16, 3);
		    $conect_c2_label_asym_id{$connect_bond} =

		      $conect_ptnr2_label_asym_id;
		    $conect_c2_auth_seq_id{$connect_bond} =

		      $conect_ptnr2_auth_seq_id;
		    $conect_c2_label_atom_id{$connect_bond} =

		      substr($connect_bond, 25, 4);
		    $conect_c2_label_alt_id{$connect_bond} =

		      substr($connect_bond, 29, 1);
		    $next_bond{++$xconnect_flag} = $connect_bond;
		}
		else {
		    $conect_symm_1{$connect_bond} = $lnks1;
		    $conect_symm_2{$connect_bond} = $lnks2;
		}
		# end if (bond_count[connect_bond]+0 == 0)

		;
	    }
	    # end for (i=1; i < link_flag; ++i)

	    ;
	}
	# end if (link_flag > 1 )
	#
	#  Process hydrogen bond information
	#
	if ($hydbnd_flag > 1) {
	    for ($i = 1; $i < $hydbnd_flag; ++$i) {
		$hbda1 = $hb_atom_beg{$i};
		$hbda2 = $hb_atom_end{$i};
		$hbdah = $hb_atom_ha{$i};
		$hbds1 = $hb_symm_1{$i};
		$hbds2 = $hb_symm_2{$i};
		$hbdsh = '   .   ';
		$hbdal1 = $hb_alt_loc_beg{$i};
		$hbdal2 = $hb_alt_loc_end{$i};
		$hbdalh = $hb_alt_loc_ha{$i};
		$hbdch = $hb_chain_id_ha{$i};
		$hbdrnh = $hb_res_seq_num_ha{$i};
		$hbdnah = $hb_name_ha{$i};
		$hbdrth = $al_res_name{($hbdch . $hbdrnh)};
		$hbdcrh = 'hydrogen';
		if (($hbdch . $hbdrnh . $hbdalh) eq ($hb_chain_id_beg{$i} .	#???

		  $hb_res_seq_num_beg{$i} . $hbdal1)) {
		    $hbrth = $hb_res_name_beg{$i};
		    $hbdsh = $hbds1;
		    $hbdcr1 = 'donor';
		    $hbdcr2 = 'acceptor';
		}
		if (($hbdch . $hbdrnh . $hbdalh) eq ($hb_chain_id_end{$i} .	#???

		  $hb_res_seq_num_end{$i} . $hbdal2)) {
		    $hbdrth = $hb_res_name_end{$i};
		    $hbdsh = $hbds2;
		    $hbdcr1 = 'acceptor';
		    $hbdcr2 = 'donor';
		}
		$hbdcr1s = $hbdcr1;
		$hbdcr2s = $hbdcr2;
		$connect_bond = ($hb_res_name_beg{$i} . $hb_chain_id_beg{$i} .

		  $hb_res_seq_num_beg{$i} . $hbda1 . $hbdal1 . ' ' .

		  $hb_res_name_end{$i} . $hb_chain_id_end{$i} .

		  $hb_res_seq_num_end{$i} . $hbda2 . $hbdal2);
		if (($hb_chain_id_beg{$i} . $hb_res_seq_num_beg{$i} . $hbda1 .

		  $hbdal1) gt ($hb_chain_id_end{$i} . $hb_res_seq_num_end{$i}	#???

		  . $hbda2 . $hbdal2)) {
		    $hbdcr1 = $hbdcr2;
		    $hbdcr2 = $hbdcr1s;
		    $hbda1 = $hb_atom_end{$i};
		    $hbda2 = $hb_atom_beg{$i};
		    $hbds1 = $hb_symm_2{$i};
		    $hbds2 = $hb_symm_1{$i};
		    $hbdal1 = $hb_alt_loc_end{$i};
		    $hbdal2 = $hb_alt_loc_beg{$i};
		    $connect_bond = ($hb_res_name_end{$i} .

		      $hb_chain_id_end{$i} . $hb_res_seq_num_end{$i} . $hbda1

		      . $hbdal1 . ' ' . $hb_res_name_beg{$i} .

		      $hb_chain_id_beg{$i} . $hb_res_seq_num_beg{$i} . $hbda2

		      . $hbdal2);
		}
		# end if ((hb_chain_id_beg[i] hb_res_seq_num_beg[i] ...
		if ($bond_count{$connect_bond} + 0 == 0) {
		    if ($hbds1 eq $hdbs2) {	#???
			$warning_list{++$warning_flag} =

			  ('#=# STRUCT_CONN: No connect for hydbnd ' .

			  $connect_bond . " \n");
		    }
		    $aslist{$connect_bond} = ('|. .|');
		    ++$bond_count{$connect_bond};
		    $conect_xtype{$connect_bond} = 'hydrog';
		    $conect_symm_1{$connect_bond} = $hbds1;
		    $conect_symm_2{$connect_bond} = $hbds2;
		    $conect_ptnr1_label_asym_id = substr($connect_bond, 4, 1);
		    $conect_ptnr1_auth_seq_id = substr($connect_bond, 5, 5);
		    $conect_ptnr2_label_asym_id = substr($connect_bond, 19,

		      1);
		    $conect_ptnr2_auth_seq_id = substr($connect_bond, 20, 5);
		    $conect_role_1{$connect_bond} = $hbdcr1;
		    $conect_role_2{$connect_bond} = $hbdcr2;
		    $conect_epsn1{$connect_bond} =

		      $ent_pol_seq_num{($conect_ptnr1_label_asym_id .

		      $conect_ptnr1_auth_seq_id)};
		    $conect_epsn2{$connect_bond} =

		      $ent_pol_seq_num{($conect_ptnr2_label_asym_id .

		      $conect_ptnr2_auth_seq_id)};
		    if ($conect_epsn1{$connect_bond} eq '') {
			$conect_epsn1{$connect_bond} = '. ';
		    }
		    if ($conect_epsn2{$connect_bond} eq '') {
			$conect_epsn2{$connect_bond} = '. ';
		    }
		    $conect_c1_label_res_id{$connect_bond} =

		      substr($connect_bond, 1, 3);
		    $conect_c1_label_asym_id{$connect_bond} =

		      $conect_ptnr1_label_asym_id;
		    $conect_c1_auth_seq_id{$connect_bond} =

		      $conect_ptnr1_auth_seq_id;
		    $conect_c1_label_atom_id{$connect_bond} =

		      substr($connect_bond, 10, 4);
		    $conect_c1_label_alt_id{$connect_bond} =

		      substr($connect_bond, 14, 1);
		    $conect_c2_label_res_id{$connect_bond} =

		      substr($connect_bond, 16, 3);
		    $conect_c2_label_asym_id{$connect_bond} =

		      $conect_ptnr2_label_asym_id;
		    $conect_c2_auth_seq_id{$connect_bond} =

		      $conect_ptnr2_auth_seq_id;
		    $conect_c2_label_atom_id{$connect_bond} =

		      substr($connect_bond, 25, 4);
		    $conect_c2_label_alt_id{$connect_bond} =

		      substr($connect_bond, 29, 1);
		    $next_bond{++$xconnect_flag} = $connect_bond;
		}
		else {
		    if ($conect_xtype{$connect_bond} ne 'hydrog') {
			$warning_list{++$warning_flag} =

			  ('#=# STRUCT_CONN: Changed bond ' . $connect_bond .

			  " to hydrog\n");
		    }
		    $conect_xtype{$connect_bond} = 'hydrog';
		    $conect_symm_1{$connect_bond} = $hbds1;
		    $conect_symm_2{$connect_bond} = $hbds2;
		    $conect_role_1{$connect_bond} = $hbdcr1;
		    $conect_role_2{$connect_bond} = $hbdcr2;
		}
		# end if (bond_count[connect_bond]+0 == 0)
		if ($hbdsh ne '   .   ' || $hbdnah ne ' .  ' ||

		  $hbdrnh ne '  .  ') {
		    $hbda2 = $hb_atom_end{$i};
		    $hbdal2 = $hb_alt_loc_end{$i};
		    $hbds2 = $hb_symm_2{$i};
		    $hbds1 = $hbdsh;
		    $hbdcr1 = $hbdcrh;
		    $hbdcr2 = $hbdcr2s;
		    $connect_bond = ($hbdrth . $hbdch . $hbdrnh . $hbdah .

		      $hbdalh . ' ' . $hb_res_name_end{$i} .

		      $hb_chain_id_end{$i} . $hb_res_seq_num_end{$i} . $hbda2

		      . $hbdal2);
		    if (($hbdch . $hbdrnh . $hbdah . $hbdalh) gt	#???

		      ($hb_chain_id_end{$i} . $hb_res_seq_num_end{$i} . $hbda2

		      . $hbdal2)) {
			$hbds1 = $hbds2;
			$hbds2 = $hbdsh;
			$hbdcr1 = $hbdcr2s;
			$hbdcr2 = $hbdcrh;
			$connect_bond = ($hb_res_name_end{$i} .

			  $hb_chain_id_end{$i} . $hb_res_seq_num_end{$i} .

			  $hbda2 . $hbdal2 . ' ' . $hbdrth . $hbdch . $hbdrnh

			  . $hbdah . $hbdalh);
		    }
		    # end if ...
		    if ($bond_count{$connect_bond} + 0 == 0) {
			$aslist{$connect_bond} = ('|. .|');
			++$bond_count{$connect_bond};
			$conect_xtype{$connect_bond} = 'hydrog';
			$conect_symm_1{$connect_bond} = $hbds1;
			$conect_symm_2{$connect_bond} = $hbds2;
			$conect_ptnr1_label_asym_id = substr($connect_bond, 4,

			  1);
			$conect_ptnr1_auth_seq_id = substr($connect_bond, 5,

			  5);
			$conect_ptnr2_label_asym_id = substr($connect_bond,

			  19, 1);
			$conect_ptnr2_auth_seq_id = substr($connect_bond, 20,

			  5);
			$conect_role_1{$connect_bond} = $hbdcr1;
			$conect_role_2{$connect_bond} = $hbdcr2;
			$conect_epsn1{$connect_bond} =

			  $ent_pol_seq_num{($conect_ptnr1_label_asym_id .

			  $conect_ptnr1_auth_seq_id)};
			$conect_epsn2{$connect_bond} =

			  $ent_pol_seq_num{($conect_ptnr2_label_asym_id .

			  $conect_ptnr2_auth_seq_id)};
			if ($conect_epsn1{$connect_bond} eq '') {
			    $conect_epsn1{$connect_bond} = '. ';
			}
			if ($conect_epsn2{$connect_bond} eq '') {
			    $conect_epsn2{$connect_bond} = '. ';
			}
			$conect_c1_label_res_id{$connect_bond} =

			  substr($connect_bond, 1, 3);
			$conect_c1_label_asym_id{$connect_bond} =

			  $conect_ptnr1_label_asym_id;
			$conect_c1_auth_seq_id{$connect_bond} =

			  $conect_ptnr1_auth_seq_id;
			$conect_c1_label_atom_id{$connect_bond} =

			  substr($connect_bond, 10, 4);
			$conect_c1_label_alt_id{$connect_bond} =

			  substr($connect_bond, 14, 1);
			$conect_c2_label_res_id{$connect_bond} =

			  substr($connect_bond, 16, 3);
			$conect_c2_label_asym_id{$connect_bond} =

			  $conect_ptnr2_label_asym_id;
			$conect_c2_auth_seq_id{$connect_bond} =

			  $conect_ptnr2_auth_seq_id;
			$conect_c2_label_atom_id{$connect_bond} =

			  substr($connect_bond, 25, 4);
			$conect_c2_label_alt_id{$connect_bond} =

			  substr($connect_bond, 29, 1);
			$next_bond{++$xconnect_flag} = $connect_bond;
		    }
		    else {
			if ($conect_xtype{$connect_bond} ne 'hydrog') {
			    $warning_list{++$warning_flag} =

			      ('#=# STRUCT_CONN: Changed bond ' .

			      $connect_bond . " to hydrog\n");
			}
			$conect_xtype{$connect_bond} = 'hydrog';
			$conect_symm_1{$connect_bond} = $hbds1;
			$conect_symm_2{$connect_bond} = $hbds2;
			$conect_role_1{$connect_bond} = $hbdcr1;
			$conect_role_2{$connect_bond} = $hbdcr2;
		    }
		    # end if (bond_count[connect_bond]+0 == 0)
		    $hbda1 = $hb_atom_beg{$i};
		    $hbdal1 = $hb_alt_loc_beg{$i};
		    $hbds2 = $hb_symm_1{$i};
		    $hbds1 = $hbdsh;
		    $hbdcr1 = $hbdcrh;
		    $hbdcr2 = $hbdcr1s;
		    $connect_bond = ($hbdrth . $hbdch . $hbdrnh . $hbdah .

		      $hbdalh . ' ' . $hb_res_name_beg{$i} .

		      $hb_chain_id_beg{$i} . $hb_res_seq_num_beg{$i} . $hbda1

		      . $hbdal1);
		    if (($hbdch . $hbdrnh . $hbdah . $hbdalh) gt	#???

		      ($hb_chain_id_beg{$i} . $hb_res_seq_num_beg{$i} . $hbda1

		      . $hbdal1)) {
			$hbds1 = $hbds2;
			$hbds2 = $hbdsh;
			$hbdcr1 = $hbdcr1s;
			$hbdcr2 = $hbdcrh;
			$connect_bond = ($hb_res_name_beg{$i} .

			  $hb_chain_id_beg{$i} . $hb_res_seq_num_beg{$i} .

			  $hbda1 . $hbdal1 . ' ' . $hbdrth . $hbdch . $hbdrnh

			  . $hbdah . $hbdalh);
		    }
		    # end if ...
		    if ($bond_count{$connect_bond} + 0 == 0) {
			$aslist{$connect_bond} = ('|. .|');
			++$bond_count{$connect_bond};
			$conect_xtype{$connect_bond} = 'hydrog';
			$conect_symm_1{$connect_bond} = $hbds1;
			$conect_symm_2{$connect_bond} = $hbds2;
			$conect_ptnr1_label_asym_id = substr($connect_bond, 4,

			  1);
			$conect_ptnr1_auth_seq_id = substr($connect_bond, 5,

			  5);
			$conect_ptnr2_label_asym_id = substr($connect_bond,

			  19, 1);
			$conect_ptnr2_auth_seq_id = substr($connect_bond, 20,

			  5);
			$conect_role_1{$connect_bond} = $hbdcr1;
			$conect_role_2{$connect_bond} = $hbdcr2;
			$conect_epsn1{$connect_bond} =

			  $ent_pol_seq_num{($conect_ptnr1_label_asym_id .

			  $conect_ptnr1_auth_seq_id)};
			$conect_epsn2{$connect_bond} =

			  $ent_pol_seq_num{($conect_ptnr2_label_asym_id .

			  $conect_ptnr2_auth_seq_id)};
			if ($conect_epsn1{$connect_bond} eq '') {
			    $conect_epsn1{$connect_bond} = '. ';
			}
			if ($conect_epsn2{$connect_bond} eq '') {
			    $conect_epsn2{$connect_bond} = '. ';
			}
			$conect_c1_label_res_id{$connect_bond} =

			  substr($connect_bond, 1, 3);
			$conect_c1_label_asym_id{$connect_bond} =

			  $conect_ptnr1_label_asym_id;
			$conect_c1_auth_seq_id{$connect_bond} =

			  $conect_ptnr1_auth_seq_id;
			$conect_c1_label_atom_id{$connect_bond} =

			  substr($connect_bond, 10, 4);
			$conect_c1_label_alt_id{$connect_bond} =

			  substr($connect_bond, 14, 1);
			$conect_c2_label_res_id{$connect_bond} =

			  substr($connect_bond, 16, 3);
			$conect_c2_label_asym_id{$connect_bond} =

			  $conect_ptnr2_label_asym_id;
			$conect_c2_auth_seq_id{$connect_bond} =

			  $conect_ptnr2_auth_seq_id;
			$conect_c2_label_atom_id{$connect_bond} =

			  substr($connect_bond, 25, 4);
			$conect_c2_label_alt_id{$connect_bond} =

			  substr($connect_bond, 29, 1);
			$next_bond{++$xconnect_flag} = $connect_bond;
		    }
		    else {
			if ($conect_xtype{$connect_bond} ne 'hydrog') {
			    $warning_list{++$warning_flag} =

			      ('#=# STRUCT_CONN: Changed bond ' .

			      $connect_bond . " to hydrog\n");
			}
			$conect_xtype{$connect_bond} = 'hydrog';
			$conect_symm_1{$connect_bond} = $hbds1;
			$conect_symm_2{$connect_bond} = $hbds2;
			$conect_role_1{$connect_bond} = $hbdcr1;
			$conect_role_2{$connect_bond} = $hbdcr2;
		    }
		    # end if (bond_count[connect_bond]+0 == 0)

		    ;
		}
	    }
	    # end for (i=1; i < hydbnd_flag; ++i)

	    ;
	}
	# end if (hydbnd_flag > 1 )
	$icid = 0;
	$model_passes = 1;
	$model_range = 0;
	$model_bias = 0;
	while ($model_bias < $model_passes) {	#???
	    for ($i = 1; $i <= $xconnect_flag; ++$i) {
		$bond = $next_bond{$i};
		if ($conect_role_1{$bond} ne '.') {
		    $conect_role_1{$bond} = ($conect_role_1{$bond} .

		      "\n                                  ");
		}
		if ($conect_role_2{$bond} ne '.') {
		    $conect_role_2{$bond} = ($conect_role_2{$bond} .

		      "\n                                                               ");
		}
		{
		    $num_aslist = (@xaslist = split(/\|/, $aslist{$bond},

		      9999));
		    if ($num_aslist >= 1 && $xaslist[$num_aslist] eq '' &&

		      '|' eq ' ') {
			--$num_aslist;
		    }
		}
		for ($ii = 1; $ii <= $num_aslist; ++$ii) {
		    if ($xaslist[$ii] ne '') {
			++$icid;
			{
			    $numxx = (@c12 = split(' ', $xaslist[$ii], 9999));
			    if ($numxx >= 1 && $c12[$numxx] eq '' &&

			      ' ' eq ' ') {
				--$numxx;
			    }
			}
			printf

			  "%6s %6s %3s %1s %5s %4s %1s %s %7s %3s %1s %5s %4s %1s %s %7s\n",

			  $icid, $conect_xtype{$bond},

			  $conect_c1_label_res_id{$bond},

			  $conect_c1_label_asym_id{$bond},

			  $conect_c1_auth_seq_id{$bond},

			  $conect_c1_label_atom_id{$bond},

			  $conect_c1_label_alt_id{$bond},

			  $conect_role_1{$bond}, $conect_symm_1{$bond},

			  $conect_c2_label_res_id{$bond},

			  $conect_c2_label_asym_id{$bond},

			  $conect_c2_auth_seq_id{$bond},

			  $conect_c2_label_atom_id{$bond},

			  $conect_c2_label_alt_id{$bond},

			  $conect_role_2{$bond}, $conect_symm_2{$bond};
			if ($model_flags eq 'yes') {
			    $cur_model = '.';
			    if ($c12[1] ne '' && $c12[1] ne '.') {
				$c12[1] = $atom_point{0 + $c12[1]};
				$test_model = int(($c12[1] + $max_asn - 1) /

				  $max_asn);
				$ihi = $model_count;
				$ilo = 1;
				if ($test_model - $ihi > 0) {
				    $test_model = $ihi;
				}
				while ($ihi - $ilo > 0) {
				    if (($c12[1] - $model_asn_low{$test_model}

				      >= 0) && ($c12[1] -

				      $model_asn_high{$test_model} <= 0)) {
					$ihi = $ilo = $test_model;
				    }
				    else {
					if ($c12[1] -

					  $model_asn_low{$test_model} < 0) {
					    $ihi = $test_model - 1;
					}
					else {
					    if ($c12[1] -

					      $model_asn_high{$test_model} >

					      0) {
						$ilo = $test_model + 1;
					    }
					    else {
						$ihi = $ilo = $test_model;
					    }
					}
					if ($ihi < $ilo) {	#???
					    $ihi = $ilo;
					}
					$test_model = int(($ihi + $ilo) / 2);
				    }
				}
				$cur_model = $model{$test_model};
				if ($model_range - $test_model < 0) {
				    $model_range = $test_model;
				}
			    }
			    if ($c12[2] ne '' && $c12[2] ne '.') {
				$c12[2] = $atom_point{0 + $c12[2]};
				if ($model_bias + 0 == 0) {
				    if (($c12[2] - $model_asn_low{$test_model}

				      < 0) || ($c12[2] -

				      $model_asn_high{$test_model} > 0)) {
					$warning_list{++$warning_flag} =

					  ('#=# STRUCT_CONN: Bond between ' .

					  $atom_number{$c12[1]} . ' and ' .

					  $atom_number{$c12[2]} .

					  " crosses models \n");
				    }
				}
			    }
			    printf

			      " %7s          %6s      %7s          %6s\n",

			      $atom_number{$c12[1]} + $model_bias *

			      $max_atom_number, $conect_epsn1{$bond},

			      $atom_number{$c12[2]} + $model_bias *

			      $max_atom_number, $conect_epsn2{$bond};
			    for ($ip = 1; $ip < 3; $ip++) {
				$id_in_model = 0;
				if ($root_model_rep{$c12[$ip]} ne '') {
				    $id_in_model =

				      $atom_id_in_model{$root_model_rep{$c12[$ip]}};
				    $auth_id_in_model =

				      $atom_number{$root_model_rep{$c12[$ip]}};
				}
				if ($id_in_model + 0 == 0) {
				    $id_in_model =

				      $atom_id_by_type{($atom_name{$c12[$ip]}

				      . '|' . $residue_name{$c12[$ip]} . '|' .

				      $atom_alt_location{$c12[$ip]} . '|' .

				      $chain_id{$c12[$ip]} . '|' .

				      $residue_seq_number{$c12[$ip]})};
				    $auth_id_in_model =

				      $atom_number{$atom_point_id_in_model{$id_in_model}};
				}
				if ($id_in_model + 0 == 0) {
				    $id_in_model = '.';
				    $auth_id_in_model = '.';
				}
				printf '   %5s           %5s     ',

				  $id_in_model, $auth_id_in_model;
			    }
			    printf "    %6s \n",

			      $model{$test_model + $model_bias};
			}
			else {
			    printf

			      "   %5s          %6s        %5s          %6s\n",

			      $c12[1], $conect_epsn1{$bond}, $c12[2],

			      $conect_epsn2{$bond};
			}
		    }
		}
	    }
	    # end for (i = 1; i <= xconnect_flag; ++i)
	    if ($model_bias + 0 == 0 && $model_flags eq 'yes' &&

	      $model_range - 1 == 0) {
		$model_passes = $model_count;
		$warning_list{++$warning_flag} =

		  ("#=# STRUCT_CONN: Replicating connectivity for all models \n");
	    }
	    $model_bias++;
	}
    }
    # end if (connect_flag+0 > 0 || ssbond_flag > 1)

    {
	if ($ss_flag > 1 && $ss_flag_2 eq '1') {	#???
	    printf (("\n\n\n"));
	    printf (("####################\n"));
	    printf (("#                  #\n"));
	    printf (("# STRUCT_CONF      #\n"));
	    printf (("#                  #\n"));
	    printf (("####################\n\n\n"));

	    printf (("\nloop_ \n"));
	    printf (("_struct_conf_type.id\n"));
	    printf (("_struct_conf_type.criteria\n"));
	    printf (("_struct_conf_type.reference\n"));

	    for ($i = 1; $i < $helix_flag; ++$i) {
		++$h_class_count{$helix_class{$i}};
		if ($h_class_count{$helix_class{$i}} == 1) {
		    printf " %s      'From PDB'   . \n", $helix_class{$i};
		}
	    }
	    if ($turn_flag > 1 && $turn_flag_2 eq '1') {	#???
		for ($i = $helix_flag; $i < $helix_flag + $turn_flag - 1;

		++$i) {
		    ++$t_class_count{$turn_class{$i}};
		    if ($t_class_count{$turn_class{$i}} == 1) {
			printf " %s      'From PDB'   . \n", $turn_class{$i};
		    }
		}
		++$turn_flag_2;
	    }

	    ++$ss_flag_2;
	}

	if ($ss_flag > 1 && $ss_flag_2 eq '2') {	#???
	    printf (("\nloop_ \n"));
	    printf (("_struct_conf.id\n"));
	    printf (("_struct_conf.conf_type_id\n"));
	    printf (("_struct_conf.beg_label_comp_id\n"));
	    printf (("_struct_conf.beg_label_asym_id\n"));
	    printf (("_struct_conf.beg_auth_seq_id\n"));
	    printf (("_struct_conf.end_label_comp_id\n"));
	    printf (("_struct_conf.end_label_asym_id\n"));
	    printf (("_struct_conf.end_auth_seq_id\n"));
	    printf (("_struct_conf.details\n"));
	    printf (("_struct_conf.beg_label_seq_id\n"));
	    printf (("_struct_conf.end_label_seq_id\n"));

	    ++$ss_flag_2;
	}

	#  start with helix records
	for ($i = 1; $i < $helix_flag; ++$i) {
	    $xxx[1] = '';
	    {
		$num_x = (@xxx = split(' ', $helix_comment{$i}, 9999));
		if ($num_x >= 1 && $xxx[$num_x] eq '' && ' ' eq ' ') {
		    --$num_x;
		}
	    }
	    $helix_comment{$i} = $xxx[1];
	    for ($j = 2; $j <= $num_x; ++$j) {
		$helix_comment{$i} = ($helix_comment{$i} . ' ' . $xxx[$j]);
	    }
	    printf "helix_%-3s %-12s %3s %1s %5s %3s %1s %5s '%s'\n",

	      $helix_id{$i}, $helix_class{$i}, $helix_res_name_beg{$i},

	      $helix_chain_id_beg{$i}, $helix_res_seq_beg{$i},

	      $helix_res_name_end{$i}, $helix_chain_id_end{$i},

	      $helix_res_seq_end{$i}, $helix_comment{$i};
	    $epsn1 = $ent_pol_seq_num{($helix_chain_id_beg{$i} .

	      $helix_res_seq_beg{$i})};
	    if ($epsn1 eq '') {	#???
		$epsn1 = '. ';
	    }
	    $epsn2 = $ent_pol_seq_num{($helix_chain_id_end{$i} .

	      $helix_res_seq_end{$i})};
	    if ($epsn2 eq '') {	#???
		$epsn2 = '. ';
	    }
	    printf "                           %6s      %6s\n", $epsn1,

	      $epsn2;
	}

	$j = $helix_flag;

	#  and now turn records
	$k = $j + $turn_flag - 1;
	for ($i = $j; $i < $k; ++$i) {
	    $xxx[1] = '';
	    {
		$num_x = (@xxx = split(' ', $turn_comment{$i}, 9999));
		if ($num_x >= 1 && $xxx[$num_x] eq '' && ' ' eq ' ') {
		    --$num_x;
		}
	    }
	    $turn_comment{$i} = $xxx[1];
	    for ($l = 2; $l <= $num_x; ++$l) {
		$turn_comment{$i} = ($turn_comment{$i} . ' ' . $xxx[$l]);
	    }
	    printf "turn_%-3s  %6s       %3s %1s %5s %3s %1s %5s '%s'\n",

	      $turn_id{$i}, $turn_class{$i}, $turn_res_name_beg{$i},

	      $turn_chain_id_beg{$i}, $turn_res_seq_beg{$i},

	      $turn_res_name_end{$i}, $turn_chain_id_end{$i},

	      $turn_res_seq_end{$i}, $turn_comment{$i};
	    $epsn1 = $ent_pol_seq_num{($turn_chain_id_beg{$i} .

	      $turn_res_seq_beg{$i})};
	    if ($epsn1 eq '') {	#???
		$epsn1 = '. ';
	    }
	    $epsn2 = $ent_pol_seq_num{($turn_chain_id_end{$i} .

	      $turn_res_seq_end{$i})};
	    if ($epsn2 eq '') {	#???
		$epsn2 = '. ';
	    }
	    printf "                           %6s      %6s\n", $epsn1,

	      $epsn2;
	}
    }
    {
	if ($cispep_flag > 0) {
	    printf (("\n\n\n"));
	    printf (("####################\n"));
	    printf (("#                  #\n"));
	    printf (("# STRUCT_MON_PROT  #\n"));
	    printf (("#                  #\n"));
	    printf (("####################\n\n\n"));

	    printf (("\nloop_ \n"));
	    printf (("_struct_mon_prot.label_comp_id\n"));
	    printf (("_struct_mon_prot.label_asym_id\n"));
	    printf (("_struct_mon_prot.auth_seq_id\n"));
	    printf (("_struct_mon_prot.label_alt_id\n"));
	    printf (("_struct_mon_prot.label_seq_id\n"));
	    printf ((('_struct_mon_prot.' . $pdb2cif_prefix .

	      "label_model_id\n")));
	    printf (("_struct_mon_prot.omega\n\n"));
	    for ($i = 1; $i <= $cispep_flag; ++$i) {
		$epsn = $ent_pol_seq_num{($cp_chain_id_end{$i} .

		  $cp_res_seq_num_end{$i})};
		if ($epsn eq '') {
		    $epsn = '. ';
		}
		printf "%s %s %s %s %s %s %6s\n", $cp_res_name_end{$i},

		  $cp_chain_id_end{$i}, $cp_res_seq_num_end{$i}, '.', $epsn,

		  $cp_modnum{$i}, $cp_omega{$i};
	    }
	    printf (("\nloop_ \n"));
	    printf (("_struct_mon_prot_cis.label_comp_id\n"));
	    printf (("_struct_mon_prot_cis.label_asym_id\n"));
	    printf (("_struct_mon_prot_cis.auth_seq_id\n"));
	    printf (("_struct_mon_prot_cis.label_alt_id\n"));
	    printf (("_struct_mon_prot_cis.label_seq_id\n"));
	    printf ((('_struct_mon_prot_cis.' . $pdb2cif_prefix .

	      "label_model_id\n")));
	    for ($i = 1; $i <= $cispep_flag; ++$i) {
		$epsn = $ent_pol_seq_num{($cp_chain_id_end{$i} .

		  $cp_res_seq_num_end{$i})};
		if ($epsn eq '') {
		    $epsn = '. ';
		}
		printf "%s %s %s %s %s %s\n", $cp_res_name_end{$i},

		  $cp_chain_id_end{$i}, $cp_res_seq_num_end{$i}, '.', $epsn,

		  $cp_modnum{$i};
	    }
	}
    }
    #
    #  Now output site information
    #
    {
	if ($site_flag_1 eq '1') {	#???
	    printf (("\n\n\n"));
	    printf (("####################\n"));
	    printf (("#                  #\n"));
	    printf (("# STRUCT_SITE      #\n"));
	    printf (("#                  #\n"));
	    printf (("####################\n\n\n"));

	    printf (("\nloop_\n"));
	    printf (("_struct_site.id\n"));
	    printf (("_struct_site.details\n"));
	    $site_flag_1 = 2;

	    for ($i = 1; $i < $site_flag; ++$i) {
		if ($site_seq_no{$i} + 0 == 1) {
		    printf " %3s  ?\n", $site_id{$i};
		}
	    }
	}

	if ($site_flag_1 eq '2') {	#???
	    printf (("\nloop_\n"));
	    printf (("_struct_site_gen.id\n"));
	    printf (("_struct_site_gen.site_id\n"));
	    printf (("_struct_site_gen.label_comp_id\n"));
	    printf (("_struct_site_gen.label_asym_id\n"));
	    printf (("_struct_site_gen.auth_seq_id\n"));
	    printf (("_struct_site_gen.label_seq_id\n"));
	    printf (("_struct_site_gen.label_alt_id\n"));
	    printf (("_struct_site_gen.symmetry\n"));
	    printf (("_struct_site_gen.details\n"));

	    $site_flag_1 = 3;
	}

	$site_num = 1;

	for ($i = 1; $i < $site_flag; ++$i) {
	    $epsn1 = $ent_pol_seq_num{($site_chain_id_1{$i} .

	      $site_res_seq_1{$i})};
	    if ($epsn1 eq '') {	#???
		$epsn1 = '.  ';
	    }
	    printf "%3s %3s %3s %1s %5s %6s . 1_555 . \n", $site_num,

	      $site_id{$i}, $site_res_name_1{$i}, $site_chain_id_1{$i},

	      $site_res_seq_1{$i}, $epsn1;
	    ++$site_num;
	    if ($site_res_name_2{$i} ne '   ') {
		$epsn2 = $ent_pol_seq_num{($site_chain_id_2{$i} .

		  $site_res_seq_2{$i})};
		if ($epsn2 eq '') {	#???
		    $epsn2 = '.  ';
		}
		printf "%3s %3s %3s %1s %5s %6s . 1_555 . \n", $site_num,

		  $site_id{$i}, $site_res_name_2{$i}, $site_chain_id_2{$i},

		  $site_res_seq_2{$i}, $epsn2;
		++$site_num;
	    }
	    if ($site_res_name_3{$i} ne '   ') {
		$epsn3 = $ent_pol_seq_num{($site_chain_id_3{$i} .

		  $site_res_seq_3{$i})};
		if ($epsn3 eq '') {
		    $epsn3 = '.  ';
		}
		printf "%3s %3s %3s %1s %5s %6s . 1_555 . \n", $site_num,

		  $site_id{$i}, $site_res_name_3{$i}, $site_chain_id_3{$i},

		  $site_res_seq_3{$i}, $epsn3;
		++$site_num;
	    }
	    if ($site_res_name_4{$i} ne '   ') {
		$epsn4 = $ent_pol_seq_num{($site_chain_id_4{$i} .

		  $site_res_seq_4{$i})};
		if ($epsn4 eq '') {
		    $epsn4 = '.  ';
		}
		printf "%3s %3s %3s %1s %5s %6s . 1_555 . \n", $site_num,

		  $site_id{$i}, $site_res_name_4{$i}, $site_chain_id_4{$i},

		  $site_res_seq_4{$i}, $epsn4;
		++$site_num;
	    }
	}
    }

    #
    #   Process SHEET information.  The PDB format makes distinct sheets
    #   out of bifurcated sheets and from sheets with broken strands.  The
    #   mmCIF format allows them to be combined.  We first must identify
    #   sheets that share strands.
    #
    #   Note:  If you convert this code to another language, be aware that
    #   strong use has been made of the awk initialization to ""
    #
    if ($sheet_flag > 0) {
	printf (("\n\n\n"));
	printf (("################\n"));
	printf (("#              #\n"));
	printf (("# STRUCT_SHEET #\n"));
	printf (("#              #\n"));
	printf (("################\n\n"));
	$x_prior_strand = '';
	$num_strand = 0;
	$num_pair = 0;
	for ($i = 1; $i <= $sheet_flag; ++$i) {
	    $strand = ($sheet_res_name_beg{$i} . ' ' . $sheet_chain_id_beg{$i}

	      . ' ' . $sheet_res_seq_beg{$i} . ' ' . $sheet_res_name_end{$i} .

	      ' ' . $sheet_chain_id_end{$i} . ' ' . $sheet_res_seq_end{$i});
	    ++$strand_count{$strand};
	    $epsn1 = $ent_pol_seq_num{($sheet_chain_id_beg{$i} .

	      $sheet_res_seq_beg{$i})};
	    $epsn2 = $ent_pol_seq_num{($sheet_chain_id_end{$i} .

	      $sheet_res_seq_end{$i})};
	    if ($epsn1 eq '') {	#???
		$epsn1 = '. ';
	    }
	    if ($epsn2 eq '') {	#???
		$epsn2 = '. ';
	    }
	    $strand_epsn1{$strand} = $epsn1;
	    $strand_epsn2{$strand} = $epsn2;
	    if ($strand_count{$strand} == 1) {
		$strand_list{++$num_strand} = $strand;
	    }
	    $my_strand_name = ($sheet_id{$i} . '|' . $sheet_strand_no{$i});
	    if ($ordered_strand_name{$strand} eq '' ||

	      $ordered_strand_name{$strand} gt $my_strand_name) {	#???
		$ordered_strand_name{$strand} = $my_strand_name;
	    }
	    if ($sheet_id_list{$strand} eq '') {
		$sheet_id_list{$strand} = $sheet_id{$i};
	    }
	    else {
		$sheet_id_list{$strand} = ($sheet_id_list{$strand} . '|' .

		  $sheet_id{$i});
	    }
	    if ($sheet_strand_no{$i} > 1) {
		$my_pair = ($strand . '|' . $x_prior_strand);
		$strand_pair{$my_pair} = $i;
		++$pair_count{$my_pair};
		if ($pair_count{$my_pair} == 1) {
		    $pair_list{++$num_pair} = $my_pair;
		}
	    }
	    $x_prior_strand = $strand;
	    if ($sheet_strands{$sheet_id{$i}} eq '') {
		++$num_sheet;
		$sheet_pdb_size{$sheet_id{$i}} = $sheet_no_strands{$i};
		$sheet_strands{$sheet_id{$i}} = $sheet_strand_no{$i};
	    }
	    else {
		$sheet_strands{$sheet_id{$i}} = ($sheet_strands{$sheet_id{$i}}

		  . ' ' . $sheet_strand_no{$i});
	    }
	}
	# prepare pointers for merge lists of sheet
	foreach $sheet_name (keys %sheet_strands) {
	    $sheet_merge{$sheet_name} = '';
	    $sheet_root{$sheet_name} = $sheet_name;
	    $sheet_size{$sheet_name} = $sheet_pdb_size{$sheet_name};
	}
	# merge sheets which share any strands
	foreach $strand (keys %sheet_id_list) {
	    {
		$num_sheet = (@sheets = split(/\|/, $sheet_id_list{$strand},

		  9999));
		if ($num_sheet >= 1 && $sheets[$num_sheet] eq '' &&

		  '|' eq ' ') {
		    --$num_sheet;
		}
	    }
	    $first_sheet = $sheets[1];
	    $first_sheet = $sheet_root{$first_sheet};
	    for ($i = 2; $i <= $num_sheet; ++$i) {
		$target_sheet = $sheets[$i];
		$target_sheet = $sheet_root{$target_sheet};
		if ($first_sheet ne $target_sheet) {	#???
		    if ($target_sheet lt $first_sheet) {	#???
			$temp_sheet = $target_sheet;
			$target_sheet = $first_sheet;
			$first_sheet = $temp_sheet;
		    }
		    $prev_sheet = $first_sheet;
		    while ($target_sheet ne '') {
			$next_sheet = $sheet_merge{$prev_sheet};
			$prev_sheet = $next_sheet;
			if ($next_sheet gt $target_sheet ||	#???

			  $next_sheet eq '') {
			    $prev_sheet = $target_sheet;
			    $target_sheet = $sheet_merge{$target_sheet};
			    $sheet_merge{$prev_sheet} = $next_sheet;
			    $sheet_root{$prev_sheet} = $first_sheet;
			    if ($sheet_size{$first_sheet} lt	#???

			      $sheet_size{$prev_sheet}) {
				$sheet_size{$first_sheet} =

				  $sheet_size{$prev_sheet};
			    }
			}
		    }
		}
	    }
	}
	# reorganize the strand names and prepare a sorted list
	$strand_point{' '} = '';
	for ($is = 1; $is <= $num_strand; ++$is) {
	    $strand = $strand_list{$is};
	    {
		$num_x = (@xxx = split(/\|/, $ordered_strand_name{$strand},

		  9999));
		if ($num_x >= 1 && $xxx[$num_x] eq '' && '|' eq ' ') {
		    --$num_x;
		}
	    }
	    $ordered_strand_name{$strand} = ($sheet_root{$xxx[1]} . ' ' .

	      $xxx[2] . ' ' . $xxx[1]);
	    {
		$numy = (@yyy = split(' ', ($xxx[2] . ' ' . $xxx[1]), 9999));
		if ($numy >= 1 && $yyy[$numy] eq '' && ' ' eq ' ') {
		    --$numy;
		}
	    }
	    $strand_name{$strand} = ($yyy[1] . '_' . $yyy[2]);
	    $prior_strand = ' ';
	    $next_strand = $strand_point{' '};
	    $reloop = 'yes';
	    while ($reloop eq 'yes') {
		if ($ordered_strand_name{$strand} lt	#???

		  $ordered_strand_name{$next_strand} || $next_strand eq '') {
		    $strand_point{$strand} = $next_strand;
		    $strand_point{$prior_strand} = $strand;
		    $reloop = 'no';
		}
		else {
		    $prior_strand = $next_strand;
		    $next_strand = $strand_point{$prior_strand};
		}
	    }
	}
	$next_strand = ' ';
	for ($i = 1; $i <= $num_strand; ++$i) {
	    $next_strand = $strand_point{$next_strand};
	    $strand_list{$i} = $next_strand;
	}
	# now sort pair names in the same order
	$pair_point{' '} = '';
	for ($ip = 1; $ip <= $num_pair; ++$ip) {
	    $my_pair = $pair_list{$ip};
	    {
		$num_x = (@xxx = split(/\|/, $my_pair, 9999));
		if ($num_x >= 1 && $xxx[$num_x] eq '' && '|' eq ' ') {
		    --$num_x;
		}
	    }
	    $x_my_pair = ($ordered_strand_name{$xxx[2]} . ' ' .

	      $ordered_strand_name{$xxx[1]});
	    $prior_pair = ' ';
	    $next_pair = $pair_point{' '};
	    $reloop = 'yes';
	    while ($reloop eq 'yes') {
		{
		    $num_y = (@yyy = split(/\|/, $next_pair, 9999));
		    if ($num_y >= 1 && $yyy[$num_y] eq '' && '|' eq ' ') {
			--$num_y;
		    }
		}
		$y_my_pair = ($ordered_strand_name{$yyy[2]} . ' ' .

		  $ordered_strand_name{$yyy[1]});
		if ($x_my_pair lt $y_my_pair || $next_pair eq '') {	#???
		    $pair_point{$my_pair} = $next_pair;
		    $pair_point{$prior_pair} = $my_pair;
		    $reloop = 'no';
		}
		else {
		    $prior_pair = $next_pair;
		    $next_pair = $pair_point{$prior_pair};
		}
	    }
	}
	$next_pair = ' ';
	for ($i = 1; $i <= $num_pair; ++$i) {
	    $next_pair = $pair_point{$next_pair};
	    $pair_list{$i} = $next_pair;
	}
	# and sort the sheet names
	$sheet_point{' '} = '';
	$num_sheet = 0;
	for ($is = 1; $is <= $num_strand; ++$is) {
	    $strand = $strand_list{$is};
	    {
		$my_num = (@my_sid = split(/\|/, $sheet_id_list{$strand},

		  9999));
		if ($my_num >= 1 && $my_sid[$my_num] eq '' && '|' eq ' ') {
		    --$my_num;
		}
	    }
	    $i = $my_sid[1];
	    $my_sheet_name = $sheet_root{$i};
	    ++$sheet_strand_count{$my_sheet_name};
	    if ($sheet_strand_count{$my_sheet_name} == 1) {
		$prior_sheet = ' ';
		$next_sheet = $sheet_point{' '};
		$reloop = 'yes';
		while ($reloop eq 'yes') {
		    if ($my_sheet_name eq $next_sheet) {	#???
			$reloop = 'no';
		    }
		    else {
			if ($my_sheet_name lt $next_sheet ||	#???

			  $next_sheet eq '') {
			    $sheet_point{$my_sheet_name} = $next_sheet;
			    $sheet_point{$prior_sheet} = $my_sheet_name;
			    $reloop = 'no';
			    ++$num_sheet;
			}
			else {
			    $prior_sheet = $next_sheet;
			    $next_sheet = $sheet_point{$prior_sheet};
			}
		    }
		}
	    }
	}
	# At this point, the sheets indexed by names are linked via
	# sheet_merge, with the root sheets given by the entries for
	# which sheet_root points to the same name and sorted lists are
	# available for sheets, strands and strand pairs
	#
	printf (("loop_\n"));
	printf (("_struct_sheet.id\n"));
	printf (("_struct_sheet.number_strands\n"));
	$my_sheet_name = ' ';
	for ($i = 1; $i <= $num_sheet; ++$i) {
	    $my_sheet_name = $sheet_point{$my_sheet_name};
	    printf "%3s %5d\n", $my_sheet_name, $sheet_size{$my_sheet_name};
	}
	printf (("\nloop_\n"));
	printf (("_struct_sheet_hbond.sheet_id\n"));
	printf (("_struct_sheet_hbond.range_id_1\n"));
	printf (("_struct_sheet_hbond.range_id_2\n"));
	printf (("_struct_sheet_hbond.range_1_beg_auth_seq_id\n"));
	printf (("_struct_sheet_hbond.range_1_beg_label_atom_id\n"));
	printf (("_struct_sheet_hbond.range_2_beg_auth_seq_id\n"));
	printf (("_struct_sheet_hbond.range_2_beg_label_atom_id\n"));
	printf (("_struct_sheet_hbond.range_1_end_auth_seq_id\n"));
	printf (("_struct_sheet_hbond.range_1_end_label_atom_id\n"));
	printf (("_struct_sheet_hbond.range_2_end_auth_seq_id\n"));
	printf (("_struct_sheet_hbond.range_2_end_label_atom_id\n"));
	printf (("_struct_sheet_hbond.range_1_beg_label_seq_id\n"));
	printf (("_struct_sheet_hbond.range_2_beg_label_seq_id\n"));
	printf (("_struct_sheet_hbond.range_1_end_label_seq_id\n"));
	printf (("_struct_sheet_hbond.range_2_end_label_seq_id\n"));
	for ($ip = 1; $ip <= $num_pair; ++$ip) {
	    $my_pair = $pair_list{$ip};
	    {
		$my_p_num = (@my_strands = split(/\|/, $my_pair, 9999));
		if ($my_p_num >= 1 && $my_strands[$my_p_num] eq '' &&

		  '|' eq ' ') {
		    --$my_p_num;
		}
	    }
	    $strand = $my_strands[1];
	    $i = $strand_pair{$my_pair};
	    $my_strand_name = $strand_name{$strand};
	    if ($sheet_strand_no{$i} > 1) {
		$p_strand = $my_strands[2];
		$p_strand_name = $strand_name{$p_strand};
		$epsn2 = $ent_pol_seq_num{($sheet_chain_id_reg_2{$i} .

		  $sheet_res_seq_reg_2{$i})};
		$epsn1 = $ent_pol_seq_num{($sheet_chain_id_reg_1{$i} .

		  $sheet_res_seq_reg_1{$i})};
		if ($epsn1 + 0 == 0) {
		    $epsn1 = '. ';
		}
		if ($epsn2 + 0 == 0) {
		    $epsn2 = '. ';
		}
		printf (("%3s %7s %7s %5s %4s %5s %4s %5s %4s %5s %4s\n" .

		  "                  %6s     %6s     %6s     %6s\n"),

		  $sheet_root{$sheet_id{$i}}, $p_strand_name, $my_strand_name,

		  $sheet_res_seq_reg_2{$i}, $sheet_atom_name_reg_2{$i},

		  $sheet_res_seq_reg_1{$i}, $sheet_atom_name_reg_1{$i},

		  $sheet_res_seq_reg_2{$i}, $sheet_atom_name_reg_2{$i},

		  $sheet_res_seq_reg_1{$i}, $sheet_atom_name_reg_1{$i}, 
		$epsn2, $epsn1, $epsn2, $epsn1);
	    }
	}
	printf (("\nloop_\n"));
	printf (("_struct_sheet_order.sheet_id\n"));
	printf (("_struct_sheet_order.range_id_1\n"));
	printf (("_struct_sheet_order.range_id_2\n"));
	printf (("_struct_sheet_order.offset\n"));
	printf (("_struct_sheet_order.sense\n"));
	for ($ip = 1; $ip <= $num_pair; ++$ip) {
	    $my_pair = $pair_list{$ip};
	    {
		$my_p_num = (@my_strands = split(/\|/, $my_pair, 9999));
		if ($my_p_num >= 1 && $my_strands[$my_p_num] eq '' &&

		  '|' eq ' ') {
		    --$my_p_num;
		}
	    }
	    $strand = $my_strands[1];
	    $i = $strand_pair{$my_pair};
	    $my_sheet_name = $sheet_root{$sheet_id{$i}};
	    $my_strand_name = $strand_name{$strand};
	    if ($sheet_strand_no{$i} > 1) {
		$p_strand = $my_strands[2];
		$j = $i - 1;
		$p_sheet_name = $sheet_root{$sheet_id{$j}};
		$p_strand_name = $strand_name{$p_strand};
		printf "%3s %7s %7s +1 %s\n", $sheet_root{$sheet_id{$i}},

		  $p_strand_name, $my_strand_name, $sheet_sense{$i};
	    }
	}

	printf (("\nloop_\n"));
	printf (("_struct_sheet_range.sheet_id\n"));
	printf (("_struct_sheet_range.id\n"));
	printf (("_struct_sheet_range.beg_label_comp_id\n"));
	printf (("_struct_sheet_range.beg_label_asym_id\n"));
	printf (("_struct_sheet_range.beg_auth_seq_id\n"));
	printf (("_struct_sheet_range.end_label_comp_id\n"));
	printf (("_struct_sheet_range.end_label_asym_id\n"));
	printf (("_struct_sheet_range.end_auth_seq_id\n"));
	printf (("_struct_sheet_range.beg_label_seq_id\n"));
	printf (("_struct_sheet_range.end_label_seq_id\n"));
	for ($is = 1; $is <= $num_strand; ++$is) {
	    $strand = $strand_list{$is};
	    {
		$my_stand_num = (@my_strand = split(' ', $strand, 9999));
		if ($my_stand_num >= 1 && $my_strand[$my_stand_num] eq '' &&

		  ' ' eq ' ') {
		    --$my_stand_num;
		}
	    }
	    for ($k = $my_strand[3]; $k <= $my_strand[6]; $k++) {	#???
		$strand_atl{($my_strand[2] . ' ' . $k)} =

		  $strand_name{$strand};
	    }
	    {
		$my_num = (@my_sid = split(/\|/, $sheet_id_list{$strand},

		  9999));
		if ($my_num >= 1 && $my_sid[$my_num] eq '' && '|' eq ' ') {
		    --$my_num;
		}
	    }
	    $i = $my_sid[1];
	    $my_sheet_name = $sheet_root{$i};
	    $my_strand_name = $strand_name{$strand};
	    $epsn1 = $strand_epsn1{$strand};
	    $epsn2 = $strand_epsn2{$strand};
	    printf "%3s %7s %s\n                %6s      %6s\n",

	      $my_sheet_name, $my_strand_name, $strand, $epsn1, $epsn2;
	}
    }

    #
    #  Prepare information for PDBX_POLY_SEQ_SCHEME
    #

    $cur_model_count = 0;
    $cur_model = '.';

    for ($i = 1; $i < $atom_flag; $i++) {
	while ($cur_model_count < $model_count &&	#???

	  $i > $model_asn_high{$cur_model_count}) {	#???
	    $cur_model_count++;
	    $cur_model = $model{$cur_model_count};
	}
	if ($cur_model_count > 1) {
	    last;
	}
	if ($atom_ent_seq_num{$i} ne '.') {
	    $seq_key = ($residue_name{$i} . '|' . $chain_id{$i} . '|' .

	      $residue_seq_number{$i} . '|' . $cur_model);
	    if ($poly_seq_al{($chain_id{$i} . ' ' . $atom_ent_seq_num{$i})} eq

	      '') {
		$poly_seq_al{($chain_id{$i} . ' ' . $atom_ent_seq_num{$i})} =

		  $i;
		$poly_seq_al_key{($chain_id{$i} . ' ' .

		  $atom_ent_seq_num{$i})} = $seq_key;
	    }
	    else {
		if ($poly_seq_al_key{($chain_id{$i} . ' ' .

		  $atom_ent_seq_num{$i})} ne $seq_key) {	#???
		    $poly_seq_al{($chain_id{$i} . ' ' .

		      $atom_ent_seq_num{$i})} = ($poly_seq_al{($chain_id{$i} .

		      ' ' . $atom_ent_seq_num{$i})} . '|' . $i);
		    $poly_seq_al_key{($chain_id{$i} . ' ' .

		      $atom_ent_seq_num{$i})} = $seq_key;
		}
	    }
	}
    }

    if ($num_poly_ents > 0) {
	printf (("\n"));
	printf (("##########################\n"));
	printf (("#                        #\n"));
	printf (("# PDBX_POLY_SEQ_SCHEME   #\n"));
	printf (("#                        #\n"));
	printf (("##########################\n\n"));
	printf (("loop_ \n"));
	printf (("_pdbx_poly_seq_scheme.asym_id\n"));
	printf (("_pdbx_poly_seq_scheme.entity_id\n"));
	printf (("_pdbx_poly_seq_scheme.seq_id\n"));
	printf (("_pdbx_poly_seq_scheme.mon_id\n"));
	printf (("_pdbx_poly_seq_scheme.auth_num\n"));
	printf (("_pdbx_poly_seq_scheme.pdb_strand_id\n"));
	for ($i = 1; $i <= $num_poly_ents; ++$i) {
	    for ($seq_number = $seq_start_num{$entity_seq_num{$entities{$i}}};

	    $seq_number <= $seq_end_num{$entity_seq_num{$entities{$i}}};	#???

	    $seq_number++) {
		$iatdat = $poly_seq_al{($entities{$i} . ' ' . $seq_number)};
		$strand_info = '.';
		if ($strand_atl{($entities{$i} . ' ' . $seq_number)} ne '') {
		    $strand_info = $strand_atl{($entities{$i} . ' ' .

		      $seq_number)};
		}
		if ($iatdat eq '') {
		    printf " %1s %5d %4s %3s    .  %s\n", $entities{$i},

		      $entity_seq_num{$entities{$i}}, $seq_number,

		      $entity_poly_seq{$seq_number}, $strand_info;
		}
		else {
		    {
			$numiatdat = (@iatdats = split(/\|/, $iatdat, 9999));
			if ($numiatdat >= 1 && $iatdats[$numiatdat] eq '' &&

			  '|' eq ' ') {
			    --$numiatdat;
			}
		    }
		    for ($k = 1; $k <= $numiatdat; $k++) {
			printf " %1s %5d %4s %3s %5s %s\n", $entities{$i},

			  $entity_seq_num{$entities{$i}}, $seq_number,

			  $entity_poly_seq{$seq_number},

			  $residue_seq_number{$iatdats[$k]}, $strand_info;
		    }
		}
	    }
	}
    }

    #
    #  print out any non-standard tokens used
    #
    if ($model_flags eq 'yes' || $seqadv_flag + 0 > 0 ||

      $connect_flag + 0 != 0 || $cispep_flag + 0 > 0) {
	printf (("\n\n\n"));
	printf (("########################\n"));
	printf (("#                      #\n"));
	printf (("# PUBL_MANUSCRIPT_INCL #\n"));
	printf (("#                      #\n"));
	printf (("########################\n\n\n"));

	printf (("\nloop_\n"));
	printf (("_publ_manuscript_incl.entry_id\n"));
	printf (("_publ_manuscript_incl.extra_item\n"));
	printf (("_publ_manuscript_incl.extra_info\n"));
	printf (("_publ_manuscript_incl.extra_defn\n\n"));

	if ($model_flags eq 'yes') {
	    $warning_list{++$warning_flag} =

	      ("#=# PUBL_MANUSCRIPT_INCL: Token _atom_site.\n" .

	      '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

	      "label_model_id used\n");
	    printf "   %s\n", $head_PDB_code;
	    printf ((("  '_atom_site." . $pdb2cif_prefix .

	      "label_model_id'\n")));
	    printf (("      'temporary place holder for NMR model ids'\n"));
	    printf (("      no\n\n"));

	    $warning_list{++$warning_flag} =

	      ("#=# PUBL_MANUSCRIPT_INCL: Token _atom_site.\n" .

	      '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

	      "id_in_model used\n");
	    printf "   %s\n", $head_PDB_code;
	    printf ((("  '_atom_site." . $pdb2cif_prefix .

	      "id_in_model'\n")));
	    printf (("      'temporary place holder for NMR id in model'\n"));
	    printf (("      no\n\n"));

	    $warning_list{++$warning_flag} =

	      ("#=# PUBL_MANUSCRIPT_INCL: Token _atom_site.\n" .

	      '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

	      "auth_id_in_model used\n");
	    printf "   %s\n", $head_PDB_code;
	    printf ((("  '_atom_site." . $pdb2cif_prefix .

	      "auth_id_in_model'\n")));
	    printf

	      (("      'temporary place holder for NMR author id in model'\n"));
	    printf (("      no\n\n"));
	}

	if ($connect_flag + 0 != 0 || $ssbond_flag > 1) {
	    $warning_list{++$warning_flag} =

	      ("#=# PUBL_MANUSCRIPT_INCL: Tokens _struct_conn.\n" .

	      '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

	      "ptnrn_atom_site_id used\n");
	    printf "   %s\n", $head_PDB_code;
	    printf ((("  '_struct_conn." . $pdb2cif_prefix .

	      "ptnr1_atom_site_id'\n")));
	    printf

	      (("      '_atom_site.id of partner 1 of structure connection'\n"));
	    printf (("      no\n\n"));
	    printf "   %s\n", $head_PDB_code;
	    printf ((("  '_struct_conn." . $pdb2cif_prefix .

	      "ptnr2_atom_site_id'\n")));
	    printf

	      (("      '_atom_site.id of partner 2 of structure connection'\n"));
	    printf (("      no\n\n"));
	    if ($model_flags eq 'yes') {
		$warning_list{++$warning_flag} =

		  ("#=# PUBL_MANUSCRIPT_INCL: Tokens _struct_conn.\n" .

		  '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

		  "ptnrn_atom_site_id_in_model used\n");
		$warning_list{++$warning_flag} =

		  ("#=# PUBL_MANUSCRIPT_INCL: Tokens _struct_conn.\n" .

		  '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

		  "ptnrn_auth_id_in_model used\n");
		$warning_list{++$warning_flag} =

		  ("#=# PUBL_MANUSCRIPT_INCL: Token _struct_conn.\n" .

		  '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

		  "label_model_id used\n");
		if ($verbose eq 'yes') {
		    $warning_list{++$warning_flag} =

		      ("#=# PUBL_MANUSCRIPT_INCL: Token _struct_conn.\n" .

		      '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

		      "auth_model_id used\n");
		}
		printf "   %s\n", $head_PDB_code;
		printf ((("  '_struct_conn." . $pdb2cif_prefix .

		  "ptnr1_id_in_model'\n")));
		printf (("      '_atom_site." . $pdb2cif_prefix .

		  "id_in_model of partner 1 of structure connection'\n"));
		printf (("      no\n\n"));
		printf "   %s\n", $head_PDB_code;
		printf ((("  '_struct_conn." . $pdb2cif_prefix .

		  "ptnr2_id_in_model'\n")));
		printf (("      '_atom_site." . $pdb2cif_prefix .

		  "id_in_model of partner 2 of structure connection'\n"));
		printf (("      no\n\n"));
		printf "   %s\n", $head_PDB_code;
		printf ((("  '_struct_conn." . $pdb2cif_prefix .

		  "ptnr1_auth_id_in_model'\n")));
		printf (("      '_atom_site." . $pdb2cif_prefix .

		  "auth_id_in_model of partner 1 of structure connection'\n"));
		printf (("      no\n\n"));
		printf "   %s\n", $head_PDB_code;
		printf ((("  '_struct_conn." . $pdb2cif_prefix .

		  "ptnr2_auth_id_in_model'\n")));
		printf (("      '_atom_site." . $pdb2cif_prefix .

		  "auth_id_in_model of partner 2 of structure connection'\n"));
		printf (("      no\n\n"));
		printf "   %s\n", $head_PDB_code;
		printf ((("  '_struct_conn." . $pdb2cif_prefix .

		  "label_model_id'\n")));
		printf (("      '_atom_site." . $pdb2cif_prefix .

		  "label_model_id of partners in structure connection'\n"));
		printf (("      no\n\n"));
		if ($verbose eq 'yes') {
		    printf "   %s\n", $head_PDB_code;
		    printf ((("  '_struct_conn." . $pdb2cif_prefix .

		      "auth_model_id'\n")));
		    printf (("      '_atom_site." . $pdb2cif_prefix .

		      "auth_model_id of partners in structure connection'\n"));
		    printf (("      no\n\n"));
		}
	    }
	}

	if ($cispep_flag + 0 != 0) {
	    $warning_list{++$warning_flag} =

	      ("#=# PUBL_MANUSCRIPT_INCL: Token _struct_mon_prot.\n" .

	      '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

	      "label_model_id used\n");
	    $warning_list{++$warning_flag} =

	      ("#=# PUBL_MANUSCRIPT_INCL: Token _struct_mon_prot_cis.\n" .

	      '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

	      "label_model_id used\n");
	    printf "   %s\n", $head_PDB_code;
	    printf ((("  '_struct_mon_prot." . $pdb2cif_prefix .

	      "label_model_id'\n")));
	    printf (("      'NMR model id'\n"));
	    printf (("      no\n\n"));
	    printf "   %s\n", $head_PDB_code;
	    printf ((("  '_struct_mon_prot_cis." . $pdb2cif_prefix .

	      "label_model_id'\n")));
	    printf (("      'NMR model id'\n"));
	    printf (("      no\n\n"));
	}

	if ($seqadv_flag > 0) {
	    $warning_list{++$warning_flag} =

	      ("#=# PUBL_MANUSCRIPT_INCL: Token _struct_ref_seq_dif.\n" .

	      '#=# PUBL_MANUSCRIPT_INCL:   ' . $pdb2cif_prefix .

	      "db_seq_num used\n");
	    printf "   %s\n", $head_PDB_code;
	    printf ((("  '_struct_ref_seq_dif." . $pdb2cif_prefix .

	      "db_seq_num'\n")));
	    printf (("      'database sequence number for alignment'\n"));
	    printf (("      no\n\n"));
	}
    }

    #
    #  print out summaries
    #

    {
	printf (("###############################################\n"));
	printf (("# This file was converted automatically from  #\n"));
	printf (("# PDB format to mmCIF format by the program   #\n"));
	printf "# pdb2cif version %-15s %11s #\n", $version, $version_date;
	printf (("#                    by                       #\n"));
	printf (("#   Phil Bourne, Herbert J. Bernstein and     #\n"));
	printf (("#            Frances C. Bernstein             #\n"));
	printf (("#                                             #\n"));
	printf (("# This work was supported in part by IUCr     #\n"));
	printf (("# (for HJB), US NSF, PHS, NIH, NCRR, NIGMS,   #\n"));
	printf (("# NLM and DOE (for FCB prior to 1998), and    #\n"));
	printf (("# and NSFgrant no. BIR 9310154 (for PEB)      #\n"));
	printf (("#                                             #\n"));

	{
	    $lx_tl = length($version);
	    $tx_tl = $version;
	    $version_lc = '';
	    for ($ix_tl = 1; $ix_tl <= $lx_tl; ++$ix_tl) {
		$cx_tl = substr($tx_tl, $ix_tl, 1);
		$cx_tl = substr(($lcaz . $cx_tl), index(($UCAZ . $cx_tl),

		  $cx_tl), 1);
		$version_lc = ($version_lc . $cx_tl);
	    }
	}

	if ($version_lc =~ /alpha/ || $version_lc =~ /beta/) {
	    printf (("# *************    WARNING   **************** #\n"));
	    printf (("# *                                         * #\n"));
	    printf (("# * THIS IS A TEST VERSION USE WITH CAUTION * #\n"));
	    printf (("# *                                         * #\n"));
	    printf (("# ******************************************* #\n"));
	}
	printf (("#                                             #\n"));
	printf (("# Conversion from PDB format to mmCIF is a    #\n"));
	printf (("# complex process.  This file should be       #\n"));
	printf (("# reviewed carefully before use.              #\n"));
	printf (("#                                             #\n"));
	printf (("# Even though the authors of pdb2cif have     #\n"));
	printf (("# made a good faith effort to ensure that     #\n"));
	printf (("# pdb2cif performs according to its           #\n"));
	printf (("# documentation, and we would greatly         #\n"));
	printf (("# appreciate hearing of any problems you      #\n"));
	printf (("# may encounter, the program pdb2cif and      #\n"));
	printf (("# any files created by pdb2cif are provided   #\n"));
	printf (("# **AS IS** without any warrantee as to       #\n"));
	printf (("# correctness, merchantability or fitness     #\n"));
	printf (("# for any particular or general use.          #\n"));
	printf (("#                                             #\n"));
	printf (("# THE RESPONSIBILITY FOR ANY ADVERSE          #\n"));
	printf (("# CONSEQUENCES FROM THE USE OF THE PROGRAM    #\n"));
	printf (("# PDB2CIF OR ANY FILE OR FILES CREATED BY IT  #\n"));
	printf (("# LIES SOLELY WITH THE USERS OF THE PROGRAM   #\n"));
	printf (("# AND FILE OR FILES AND NOT WITH AUTHORS OF   #\n"));
	printf (("# PDB2CIF                                     #\n"));
	printf (("#                                             #\n"));
	printf (("# The program pdb2cif is available from       #\n"));
	printf (("# http://www.bernstein-plus-sons.com          #\n"));
	printf (("#                      /software/pdb2cif or   #\n"));
	printf (("# the IUCr and its mirrors (see               #\n"));
	printf (("# http://www.iucr.org/iucr-top/cif            #\n"));
	printf (("#                      /software/pdb2cif) or  #\n"));
	printf (("# SDSC (see                                   #\n"));
	printf (("# http://www.sdsc.edu/pb/pdb2cif/pdb2cif) or  #\n"));
	printf (("# NDB and its mirrors (see                    #\n"));
	printf (("# http://ndbserver.rutgers.edu/mmcif          #\n"));
	printf (("#                      /software/pdb2cif)     #\n"));
	printf (("# and the NDB mirror sites                    #\n"));
	printf (("#                                             #\n"));

	printf (("#    See H. Bernstein, F. Bernstein,          #\n"));
	printf (("#    P. E. Bourne CIF Applications. VIII.     #\n"));
	printf (("#    pdb2cif: Translating PDB Entries into    #\n"));
	printf (("#    mmCIF Format, J. Appl. Cryst., 31,       #\n"));
	printf (("#    1998, pp 282-295.                        #\n"));
	printf (("#                                             #\n"));

	$my_at = '@';
	printf (("# Please report problems to:                  #\n"));
	printf ((('#         pdb2cif' . $my_at .

	  "bernstein-plus-sons.com     #\n")));
	printf (("#                                             #\n"));
	printf (("###############################################\n\n\n"));

	printf "\n# REMARK records parsed \t= %6d;", $all_remarks;
	printf "# specified by PDB \t= %6d\n", $total_remark;

	for ($X = 1; $X <= 999; ++$X) {
	    $total_ftnote_flag = $total_ftnote_flag + $ftnote_flag{$X};
	}
	# end for (x = 1; x <= 999; ++x)
	printf "# FTNOTE records parsed \t= %6s;", $total_ftnote_flag;
	printf "# specified by PDB \t= %6s\n", $total_ftnote;

	printf "# HET records parsed \t\t= %6d;", ($het_flag - 1);
	printf "# specified by PDB \t= %6d\n", $total_het;

	printf "# HELIX records parsed \t\t= %6d;", ($helix_flag - 1);
	printf "# specified by PDB \t= %6d\n", $total_helix;

	printf "# SHEET records parsed \t\t= %6d;", ($sheet_flag);
	printf "# specified by PDB \t= %6d\n", $total_sheet;

	printf "# TURN records parsed \t\t= %6d;", ($turn_flag - 1);
	printf "# specified by PDB \t= %6d\n", $total_turn;

	printf "# SITE records parsed \t\t= %6d;", ($site_flag - 1);
	printf "# specified by PDB \t= %6d\n", $total_site;

	$total_a_h_flag = ($atom_flag - 1 - $ter_flag);

	printf "# AT+HET records parsed \t= %6d;", ($total_a_h_flag);
	printf "# specified by PDB \t= %6d\n", $total_a_h;

	printf "# TER records parsed \t\t= %6d;", $ter_flag;
	printf "# specified by PDB \t= %6d\n", $total_ter;

	printf "# CONECT records parsed \t= %6d;", $conect_flag_2;
	printf "# specified by PDB \t= %6d\n", $total_conect;

	printf "# SEQRES records parsed \t= %6d;", ($seqres_flag - 1);
	printf "# specified by PDB \t= %6d\n", $total_seqres;
	printf "# Total of %6d records processed from PDB file\n", $.;
	for ($i = 1; $i <= $warning_flag; ++$i) {
	    printf (($warning_list{$i}));
	}
    }
}

