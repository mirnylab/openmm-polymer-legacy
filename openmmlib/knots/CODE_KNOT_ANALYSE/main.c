/**********************************
 * main.c                         *
 *                                *
 * last modification : 04/19/2004 *
 **********************************/


/* include header files */
#include "element.h"
#include "r250.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>


/* global variables */
int         nrmonA, nrchainsA=1, polA, mcs=0;
MYVEC       *posA, *cm0A;
DATA        data; /* allocation ? */

int         nrmon_knotA, nrmon_knotB;
MYVEC       *knotA, *knotB, *closure_mon;

int         distance_to_end;


static char *filenamein  = NULL;
char *redknotout = NULL;
int   redknotcnt = 0;
double knotAt = -1.1;

#define USAGE "%s [-r <reducedout>] filein\n"
#define ARGS "r:f:p:"

int parse_args(int argc, char **argv){
  int c,errflg = 0;

  optarg = NULL;
  optind = 0;

  while(!errflg && (c=getopt(argc,argv,ARGS))!=-1){
    switch(c){
      case 'r':
	redknotout = optarg;
	break;
      case 'f':
	filenamein  = optarg;
	break;
      case 'p':
	knotAt = strtod(optarg, NULL);
	printf("Alexanders_at_%f_", knotAt);
	break;
    }
  }

  if(optind < 1 || argc < 2)
    errflg=1;

  if(optind < argc)
      filenamein = argv[optind];

  if(errflg || filenamein == NULL){
    fprintf(stderr,USAGE, argv[0]);
    exit(0);
  }



  return 0;
}





int main(int argc, char **argv) {
  Py_SetProgramName (argv [0]);
  Py_Initialize ();
  initdet();


  int i=0, seed=-1;
  //char filenamein[100]="/var/www/cgi-bin/knots/tmp/conf_in", filenameout[30];
  
  parse_args(argc, argv);

  /*-----------------------------------*
   | read monomer vectors from conf_in |
   *-----------------------------------*/

  sysin(filenamein);    
  if(nrmonA!=0){
    check_for_knots(knotAt);
  }
  else {printf("0\t0\t0\t0\t0\t0\n");}

  /*---------------------------------------------*
   | initialize seed for random number generator |
   *---------------------------------------------*/

  /* if(seed<0) time ((time_t *)&seed);
  r250_srandom(seed); */

  /*---------------------------------------------*
   | generate 3d random walk and check for knots |
   *---------------------------------------------*/

  /* while(mcs<12500000){
    generate_3d_random_walk();
    //if(mcs==9213){
    check_for_knots();//}
    mcs++;
  } */

  /* sprintf(filenameout,"conf_out");
  sysout(filenameout,posA,nrmonA,mcs);  */
  Py_Finalize ();
  return(OK);
}





