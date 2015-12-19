/*************************************************************
 * analysis.c: analysis-fcts update data-array before output *
 *             file "histo.dat" is written                   *
 *                                                           * 
 * last modification: 04/19/2004                             *  
 *************************************************************/


#include "element.h"
#include <math.h>
#include <stdio.h>

extern int       nrchainsA, nrmonA, polA;
extern MYVEC     *posA, *cm0A;

extern DATA      data;


/*******************************************
 * analysis: updates data array by calling *
 *           most analysis functions,      *
 *           updates "histo.dat" output    *
 *           file                          *
 *                                         *
 * last modification: 04/19/2004           *
 *******************************************/

void analysis(int zeit){      

  FILE   *datei;
  char   fname[80];

  /*-------------------------*
   | call analysis functions |
   *-------------------------*/

  bondwinkel();        /* calculates bondangles               */
  bondlength();        /* bondlength, b**2, bmax              */
  calc_polymer();      /* end-to-end and radius of gyration   */
  calc_msd();          /* calculates mean square displacement */

  sprintf(fname,"histo.dat");

  if ((datei=fopen(fname,"a"))==NULL) {
    fprintf(stderr,"ERROR: Couldn't open %s\n",fname);
    exit(ERROR); } 

  /*------------------------------*
   | write results to "histo.dat" |
   *------------------------------*/

  fprintf(datei,"%d %lf %lf %lf %lf %lf\n",zeit,data.bondangleA,data.bondlength1A,data.endA,data.gyrA,data.msdA);

  fclose(datei); 
}


/*****************************************
 * bondwinkel: calculates cosine of      *
 *             bond-angle and averages   *
 *             over all bonds            *
 *                                       *
 * last modification: 04/19/2004         *
 *****************************************/

void bondwinkel() {
  int i,j;
  double awinkel;
  MYVEC *mon;
  MYVEC bond[2];
  
  /*----------------------------------*
   | calculate average cos(bondangle) |  
   *----------------------------------*/

  awinkel=0;

  for(i=0;i<nrchainsA;i++){
    
    mon=posA+i*polA;

    for(j=0;j<(polA-2);j++,mon++) {
      DIFF(mon+1,mon,bond);
      DIFF(mon+1,mon+2,bond+1);
      awinkel+= CWINKEL(bond,bond+1);  
      /* cos(bondangle) ! - see def. element.h */
    }
  }

  if (nrchainsA>0)
    awinkel /= ((double)(nrmonA-2*nrchainsA));
  else
    awinkel=0;
  
  /* update data array */
  data.bondangleA = awinkel;
}


/*****************************************
 * bondlength: calculates av.bondlenght, *
 *             b**2, bmax                *
 *                                       *
 * last modification: 04/19/2004         *
 *****************************************/

void bondlength() {    
  int i,j;
  double length,abslength,maxlength,distance;
  MYVEC *mon;


  /*---------------------------------------------*
   | calculate average bondlength, b**2 and bmax |
   *---------------------------------------------*/

  length=0;    /* average bondlength */
  abslength=0; /* (average bondlength)^2 */ 
  maxlength=0; /* maximum bondlength */

  for(i=0;i<nrchainsA;i++){
    
    mon=posA+i*polA;

    for(j=0;j<(polA-1);j++,mon++) {

      distance   = DIST(mon,mon+1);
      length    += (distance);
      abslength += SQUARE(distance);
      if (ABS(distance)>maxlength) maxlength=ABS(distance);

      if (ABS(distance)>2.0) {                             
      /* check if chain is broken */                          

        fprintf(stderr,"WARNING: typeA: chain %d mon "
                       "%i\nlaenge=%f\n",i,j,ABS(distance));

        fprintf(stderr,"WARNING: monomer 1   %f %f %f\n",
                       mon->x,mon->y,mon->z);
        fprintf(stderr,"WARNING: monomer 2   %f %f %f\n",
                       (mon+1)->x,(mon+1)->y,(mon+1)->z);
      }
    }
  }
  
  if ((nrchainsA>0)&&(polA>1))
    data.bondlength1A = length/((double)(nrmonA-nrchainsA));
  else
    data.bondlength1A = 0;
    
  if ((nrchainsA>0)&&(polA>1))
    data.bondlength2A = abslength/((double)(nrmonA-nrchainsA));
  else
    data.bondlength2A = 0;

  data.bondlengthmaxA = maxlength;
}


/*****************************************
 * calc_polymer: calculates end-to-end   *
 *               distance and radius of  *
 *               gyration                *
 *                                       *
 * last modification: 04/19/2004         *
 *****************************************/
  
void calc_polymer(){

  double end,gyration;
  MYVEC *mon1,*mon2;
  int i,j,k;


  /*---------------------------------------*
   | calculate end-to-end distance and     |
   | radius of gyration                    |
   *---------------------------------------*/

  end      = 0;
  gyration = 0;

  for(i=0;i<nrchainsA;i++){

    mon1 = posA + i*polA;
    end += DISTANCE2(mon1,mon1+polA-1);

    for(j=0;j<polA;j++,mon1++) {
      mon2 = posA +i*polA;
      for(k=0;k<polA;k++,mon2++)
        gyration += DISTANCE2(mon1,mon2);
    }
  }
  if (nrchainsA>0) {
    end      /= nrchainsA;
    gyration /= (2.0*SQUARE(polA)*nrchainsA);
  }
  else {
    end      = 0;
    gyration = 0;
  }
 
  data.endA = end;
  data.gyrA = gyration;
}


/******************************************
 * calc_cm: calculates center-of-mass of  *
 *          polymer at position "start"   *
 *          in pos-array                  *
 *                                        *
 * last modification: 04/19/2004          *
 ******************************************/

void calc_cm(int pol, MYVEC *start, MYVEC *result){

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


/******************************************
 * calc_msd: calculates mean-square-      *
 *           displacement                 *
 *                                        *
 * last modification: 04/19/2004          *
 ******************************************/

void calc_msd() {

  int i;
  MYVEC cm, *mon, *vhilf;
  double rmove;


  /*------------------------------------*
   | calculate mean-square-displacement |
   *------------------------------------*/ 

  vhilf= cm0A;
  rmove=0;

  for(i=0;i<nrchainsA;i++,vhilf++) {

    mon = posA + i*polA;
    calc_cm(polA,mon,&cm);

    rmove += (SQUARE((vhilf->x - cm.x))+\
              SQUARE((vhilf->y - cm.y))+\
              SQUARE((vhilf->z - cm.z)));
  }
 
  if (nrchainsA>0)
    data.msdA = rmove/nrchainsA;
  else
    data.msdA = 0;
}
