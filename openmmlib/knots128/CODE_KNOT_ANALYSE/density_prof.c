/*********************************************************
 * density_prof.c: determines density profile of globule *
 *                                                       *
 * last modification: 07/09/2004                         *
 *********************************************************/


#include "element.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern int mcs, polA;
extern int nrmonA;
extern MYVEC *posA;


/******************************************
 * calc_cm1: calculates center-of-mass of *
 *          polymer at position "start"   *
 *          in pos-array                  *
 *                                        *
 * last modification: 07/09/2004          *
 ******************************************/
                                                                                                                             
void calc_cm1(int pol, MYVEC *start, MYVEC *result){
                                                                                                                             
  int j;
  MYVEC *mon;
                                                                                                                             
  result->x = 0;
  result->y = 0;
  result->z = 0;
                                                                                                                             
  mon = start;
                                                                                                                             
  for(j=0;j<pol;j++,mon++) {
    result->x += mon->x;
    result->y += mon->y;
    result->z += mon->z;
  }
                                                                                                                             
  result->x /= pol;
  result->y /= pol;
  result->z /= pol;
}


/***************************************************
 * dens_prof: calculates density profile and       *
 *            stores info in file "density_prof"   *
 *                                                 *
 * input: *pos:    pointer to first element in     *
 *                 monomer-vector array            *
 *        nrmon:   # of monomers in that array     *
 *        abstand: distance between "Kugelschalen" *
 *                                                 *
 * last modification: 07/14/2004                   *
 ***************************************************/

void dens_prof(MYVEC *pos, int nrmon, double abstand){

  MYVEC cm;
  int spheres [1000];
  int i;
  double distance;
  MYVEC *c0,*c1;

  FILE   *datei;
  char   fname[80];


  /*----------------------------------------------------------------------*
   | 1. initialize array and calculate center of mass of complete globule |
   *----------------------------------------------------------------------*/

  for(i=0;i<1000;i++){
    spheres[i]=0;
  }
  calc_cm1(nrmonA, posA, &cm); 

  /*-----------------------------------------------------*
   | 2. put particles into spheres around center of mass |
   *-----------------------------------------------------*/

  /* for(i=0;i<nrmon;i++) { 
    distance=DIST(pos+i,&cm);
    spheres[(int)(distance/abstand)]++;
  } */

  /*--------------------------*
   | 3. print results to file |
   *--------------------------*/

  /* sprintf(fname,"density_profile");
  datei=fopen(fname,"a");
 
  for(i=0;i<100;i++){ 
    fprintf(datei,"%d\t%d\n",i+1,spheres[i]); 
  }
  fclose(datei); */

  /*-------------------------------------*
   | 4. determine "profile" of endpoints |
   *-------------------------------------*/

  sprintf(fname,"new_coeff1");
  datei=fopen(fname,"a");
                                                                                                                             
  i=(NMONOMAXA/2-1)+0.1*NMONOMAXA;
  fprintf(datei,"%d\t1\n",(int)((DIST(pos+i,&cm)/abstand)+1) );
  /* fprintf(datei,"%lf\t",DIST(pos+i,&cm) ); */

  i=(NMONOMAXA/2-1)-0.1*NMONOMAXA;
  fprintf(datei,"%d\t1\n",(int)((DIST(pos+i,&cm)/abstand)+1) ); 
  /* fprintf(datei,"%lf\t",DIST(pos+i,&cm) ); */

  /* fprintf(datei,"%lf\t",DIST(pos,pos+i) ); */

  fclose(datei);
}

