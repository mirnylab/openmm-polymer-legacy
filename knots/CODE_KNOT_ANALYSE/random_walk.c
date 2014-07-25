/****************************************************
 * random_walk.c: functions to generate and monitor *
 *                a 3d random walk                  *
 *                                                  *
 * last modification: 04/19/2004                    *
 ****************************************************/

#include "element.h"
#include "r250.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NRND_MAX (2*NMONOMAXA+1)

extern int     polA, nrmonA;
extern MYVEC   *posA;


/*****************************************************************
 * generate_3d_random_walk(): generates 3d random walk of length *
 *                            NMONOMAX from scratch and saves it *
 *                            in posA-array                      *
 *                                                               *
 * returns: OK                                                   *
 * 		                                                 * 
 * last modification: 04/19/2004                                 *
 *****************************************************************/

int generate_3d_random_walk()
{

  static double *rndtab = NULL;
  static double *rnd_ptr;
  int acc;

  int i;  
  double length;

  MYVEC *mon, *pos;
  double ctheta,stheta,phi;
  
  /*------------------------------------------*
   | initialize or update random number array |
   *------------------------------------------*/
 
  if (rndtab==NULL) {  /* first call -> initialization */ 

  /* allocate memory for NRND_MAX #s */     
  rndtab = (double*) calloc(NRND_MAX,sizeof(double));

  /* assign random #s to array       */   
  double_r250_vector(rndtab,NRND_MAX);

  /* set rnd_ptr to first position   */
  rnd_ptr = rndtab;
  }

  acc = (int)(rnd_ptr-rndtab);
    if ( (acc<0)||(acc>NRND_MAX) ) {
      fprintf(stderr,"ERROR: must generate too many RNDs"
                   "in rep.c (%d)\n",acc);
  }                                                                                                         

  else {   /* update random number array      */
  
  /* refresh rnd #s which have already been used */
  double_r250_vector(rndtab,acc);
  
  /* set rnd_ptr to first position */
  rnd_ptr = rndtab;
  }
  

  /*--------------------------------*
   | set first "particle" to origin |
   *--------------------------------*/

  mon=posA;
 
  mon->x = 0;
  mon->y = 0;
  mon->z = 0;

  nrmonA = 1;
  mon->monomer_number = nrmonA;
  mon++;
  pos=posA;

  /*------------------------------*
   | generate rest of random walk |
   *------------------------------*/

  for(i=1;i<NMONOMAXA;i++){

    length = 1;

    ctheta = 2*(*rnd_ptr++) - 1;
    stheta = 2*(*rnd_ptr++) - 1;
                                                                                                            
    phi    = ctheta*ctheta + stheta*stheta;
    while (phi >= 1 || phi == 0) {
      ctheta = 2*double_r250() - 1;
      stheta = 2*double_r250() - 1;
      phi    = ctheta*ctheta + stheta*stheta;
    }

    mon->z = pos->z + length*(1.0 - 2*phi);
    phi    = 2*sqrt(1-phi)*length;
    mon->x = pos->x + phi*ctheta;
    mon->y = pos->y + phi*stheta;
      /* if (mon->z<0 || ( (mon->z) < 0.577350269*sqrt((mon->x)*(mon->x)) ) ) return(0); */  /* wedge */

  nrmonA++;
  mon->monomer_number = nrmonA;
  mon++;
  pos++;
  }

  return(1);
}


int setup_loop()
{
  int i; 
  double s,t,x,y;

  nrmonA=0;
  s=PI/(NMONOMAXA-1);
  for(i=0;i<=(NMONOMAXA-1);i++){
    t=2*i*s;
    (posA+i)->monomer_number=i+1;
    (posA+i)->x=1/(2*sin(s))*sin(t);
    (posA+i)->y=1/(2*sin(s))*cos(t);
    (posA+i)->z=0; 
    nrmonA++;
    /* printf("%.20lf\t%.20lf\n",x,y)*/;
  }
  
  return(0);
}


/*****************************************************************
 * rotate_loop(): - selects two particles at random              *
 *                - rotates part of loop, which is confinded by  *
 *                  these particles, by random angle (0-2*Pi)    *
 *                                                               *
 * returns: OK                                                   *
 *                                                               *
 * last modification: 02/25/2004                                 *
 *****************************************************************/
 
int rotate_loop()
{
 
  static double *rndtab = NULL;
  static double *rnd_ptr;
  int acc;
 
  int i;
  int ident1,ident2,start,stop;

  MYVEC *mon, *W, *V, *V_new, *V_cross_W; 
  MYVEC mon_d, W_d, V_d, V_new_d, V_cross_W_d;
  double phi, sin_phi, cos_phi;
  double norm;

  mon=&mon_d;
  W=&W_d;
  V=&V_d;
  V_new=&V_new_d;
  V_cross_W=&V_cross_W_d;
   
  /*------------------------------------------*
   | initialize or update random number array |
   *------------------------------------------*/
  
  if (rndtab==NULL) {  /* first call -> initialization */
 
  /* allocate memory for NRND_MAX #s */
  rndtab = (double*) calloc(NRND_MAX,sizeof(double));
 
  /* assign random #s to array       */
  double_r250_vector(rndtab,NRND_MAX);
 
  /* set rnd_ptr to first position   */
  rnd_ptr = rndtab;
  }
 
  acc = (int)(rnd_ptr-rndtab);
    if ( (acc<0)||(acc>NRND_MAX) ) {
      fprintf(stderr,"ERROR: must generate too many RNDs"
                   "in rep.c (%d)\n",acc);
  }
 
  else {   /* update random number array      */
   
  /* refresh rnd #s which have already been used */
  double_r250_vector(rndtab,acc);
   
  /* set rnd_ptr to first position */
  rnd_ptr = rndtab;
  }
   
   
  /*----------------------------------------------------------------------*
   | 1. -select two particles at random around which loop will be rotated |
   |    -select rotation angle (0<=phi<2*Pi)                              |
   *----------------------------------------------------------------------*/
 
  ident1  = (int) ( ((double)(nrmonA))*(*rnd_ptr++) );
     /* choose monomer */
     if(ident1==(nrmonA)) ident1-=1;  /* if rnd_number is exactly 1 -> possible segmentation fault ! */
 
  ident2  = (int) ( ((double)(nrmonA))*(*rnd_ptr++) );
     if(ident2==(nrmonA)) ident2-=1;
 
  if(ident1==ident2) return(1);
  else if(ident1>ident2) {
      start=ident2;
      stop=ident1;
  }
  else if(ident2>ident1) {
      start=ident1;
      stop=ident2;
  }

  /* determine rotation angle */
  phi = (2*PI)*(*rnd_ptr++);
  /* printf("%lf\n",phi);*/
 
  /* W: rotation axis */
  DIFF(posA+stop,posA+start,mon);
  norm=1/BETRAG(mon);
  SDOT(norm,mon,W); /* normalize W */
  /* printf("%lf\n",BETRAG(W));*/
  cos_phi=cos(phi);
  sin_phi=sin(phi);
 
 
  /*----------------*
   | 2. rotate loop |
   *----------------*/
 
  for(i=start+1;i<stop;i++){
 
    /* V: vector to be rotated */
    DIFF(posA+i,posA+start,V);
    VDOT(V,W,V_cross_W);
 
    /* rotate loop */
    V_new->x = (V->x)*cos_phi+(1-cos_phi)*DOT(V,W)*(W->x)-(V_cross_W->x)*sin_phi;
    V_new->y = (V->y)*cos_phi+(1-cos_phi)*DOT(V,W)*(W->y)-(V_cross_W->y)*sin_phi;
    V_new->z = (V->z)*cos_phi+(1-cos_phi)*DOT(V,W)*(W->z)-(V_cross_W->z)*sin_phi;
 
    /* determine new position of monomer */
    ADD(posA+start,V_new,posA+i);
  }

  return(0);
}
