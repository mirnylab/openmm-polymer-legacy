/*********************************************************
 * knot.c: functions to determine knots in configuration *
 *                                                       *
 * last modification: 03/25/2004                         *
 *********************************************************/


#include "element.h"
#include <stdio.h>
#include <stdlib.h>
#include "r250.h"

extern int     mcs, polA, nrmonA;
extern MYVEC   *posA;
extern int     nrmon_knotA, nrmon_knotB;
extern MYVEC   *knotA, *knotB, *closure_mon;

extern DATA    data;


/************************************************************
 * check_for_knots(): - reduces complexity of configuration *
 *                      by calling reduce_complexity        *
 *                    - analyzes remaining structure with   *
 *                      knot polynomials                    *
 *                                                          *    
 * returns:	      YES: found a knot                     *
 * 		      NO : didn't find a knot               *
 *							    *		
 * last modification: 11/30/2004                            *
 ************************************************************/

int check_for_knots(double knotAt)
{
  int i,j=0,ctr=0;
  int start=1, end=nrmonA-2, size, position, new_first_position;
  MYVEC *mon, *mon_knotB;
  char fname[30];
  double t1,t2; 
  long double polynomial, polynomial1, polynomial2;
  double bondlength;
  int fixed1, fixed2;
  FILE   *datei;

  /* allocate memory for arrays which are needed to analyze knots */
  knotA           = (MYVEC *) calloc(NMONOMAXA+10, sizeof(MYVEC));
  knotB           = (MYVEC *) calloc(NMONOMAXA+10, sizeof(MYVEC)); 
  closure_mon     = (MYVEC *) calloc(2, sizeof(MYVEC));
 
  polynomial=test_for_knots (0, 0, knotAt, 1, nrmonA-2, &start, &end);
  // printf("%d\t%.5lf\t",nrmonA,polynomial); 

    printf("%d\t%.5Lf\t",nrmonA,polynomial);
  //}

  //else {printf("0\t0\t0\t0\n");}

  free(knotA);
  free(knotB);
  free(closure_mon);
  return(0);
}


/********************************************************************
 * closure: closes open polymer such that ring is formed            *
 *                                                                  *
 * int method: 0 close loop by connecting first and last particle   *
 *             1 " " " " first and last particle to infinity        *
 *               (along x-axis)                                     *
 *             2 like 1 but connection along line through endpoints *
 *                                                                  *
 * last modification: 10/12/2004                                    *
 ********************************************************************/

void closure (int method)
{
MYVEC mon[1],vhelp[1];
double large_value1, large_value2;

  if(method==0) {
    VCOPY(posA,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;
  }

  if(method==10) {
    VCOPY(knotB,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;
  }
 
  if(method==11) {
    mon->x=1.1;
    mon->y=1.1;
    mon->z=1000000.1;
    VCOPY(mon,knotB+nrmon_knotB) 
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;

    VCOPY(knotB,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+2;
    nrmon_knotB++;
  }

  if(method==1) {

    if( (knotB+nrmon_knotB-1)->x > knotB->x) {
    large_value1=100000000;
    large_value2=100000000;
    }
    else {
    large_value1=-100000000;
    large_value2=100000000;
    }

    mon->x = ((knotB+nrmon_knotB-1)->x)+large_value1;
    mon->y = ((knotB+nrmon_knotB-1)->y);
    mon->z = ((knotB+nrmon_knotB-1)->z); 
       
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;

    mon->x = ((knotB+nrmon_knotB-2)->x)+large_value1;
    mon->y = ((knotB+nrmon_knotB-2)->y)-large_value2;
    mon->z = ((knotB+nrmon_knotB-2)->z);
                                                                                                                             
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+2;
    nrmon_knotB++;

    mon->x = ((knotB)->x)-large_value1;
    mon->y = ((knotB)->y)-large_value2;
    mon->z = ((knotB)->z);
                                                                                                                             
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+3;
    nrmon_knotB++;

    mon->x = ((knotB)->x)-large_value1;
    mon->y = ((knotB)->y);
    mon->z = ((knotB)->z);
                                                                                                                             
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+4;
    nrmon_knotB++;

    mon->x = ((knotB)->x);
    mon->y = ((knotB)->y);
    mon->z = ((knotB)->z);
                                                                                                                            
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+5;
    nrmon_knotB++;
  }

  if(method==2){
  
    large_value1 = 100000000;
    large_value2 = large_value1/DIST(knotB+nrmon_knotB-1,knotB);

    //if(large_value2<1000) large_value2 = 100000000;

    //VPRINT(knotB+nrmon_knotB-1);
    //VPRINT(knotB);
    //printf("%lf\n",large_value2);

    mon->x = knotB->x + large_value2 * ((knotB+nrmon_knotB-1)->x - knotB->x);
    mon->y = knotB->y + large_value2 * ((knotB+nrmon_knotB-1)->y - knotB->y);
    mon->z = knotB->z + large_value2 * ((knotB+nrmon_knotB-1)->z - knotB->z);

    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    //VPRINT(knotB+nrmon_knotB); 
    nrmon_knotB++;

    DIFF(knotB+nrmon_knotB-2,knotB,mon); 
    VDOT(posA+1,mon,vhelp);
    //VDOT(knotB+1,mon,vhelp);
    //VPRINT(vhelp);
    mon->x = knotB->x + large_value1 * (vhelp->x)/BETRAG(vhelp);
    mon->y = knotB->y + large_value1 * (vhelp->y)/BETRAG(vhelp);
    mon->z = knotB->z + large_value1 * (vhelp->z)/BETRAG(vhelp);  
  
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+2;
    nrmon_knotB++;

    mon->x = knotB->x - large_value2 * ((knotB+nrmon_knotB-3)->x - knotB->x);
    mon->y = knotB->y - large_value2 * ((knotB+nrmon_knotB-3)->y - knotB->y);
    mon->z = knotB->z - large_value2 * ((knotB+nrmon_knotB-3)->z - knotB->z);

    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+3;
    nrmon_knotB++;
 
    VCOPY(knotB,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+4;
    nrmon_knotB++;
  }
 
  if(method==3){

    calc_cm(nrmon_knotB, knotB, vhelp);  /* save center of mass in vhelp */
    large_value1 = 100000000;
  
    large_value2 = large_value1/DIST(vhelp,knotB+nrmon_knotB-1);
    mon->x = vhelp->x + large_value2 * ((knotB+nrmon_knotB-1)->x - vhelp->x);
    mon->y = vhelp->y + large_value2 * ((knotB+nrmon_knotB-1)->y - vhelp->y);
    mon->z = vhelp->z + large_value2 * ((knotB+nrmon_knotB-1)->z - vhelp->z); 

    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;

    large_value2 = large_value1/DIST(vhelp,knotB);
    mon->x = vhelp->x + large_value2 * (knotB->x - vhelp->x);
    mon->y = vhelp->y + large_value2 * (knotB->y - vhelp->y);
    mon->z = vhelp->z + large_value2 * (knotB->z - vhelp->z);
 
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+2;
    nrmon_knotB++;

    VCOPY(knotB,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+3;
    nrmon_knotB++;
  }
}


void change_posA()
{
                                                                                                                             
int i,j;
double large_value1,large_value2;
MYVEC mon1[1],mon2[1],vhelp[1];
 
    copy_MYVEC_array(posA,knotB,0,nrmonA);
 
    calc_cm(nrmonA, posA, vhelp);  /* save center of mass in vhelp */
    large_value1 = 100000000;
 
    // calculate first
    large_value2 = large_value1/DIST(vhelp,posA);
    mon1->x = vhelp->x + large_value2 * (posA->x - vhelp->x);
    mon1->y = vhelp->y + large_value2 * (posA->y - vhelp->y);
    mon1->z = vhelp->z + large_value2 * (posA->z - vhelp->z);                                                                                                                              
    // calculate last
    large_value2 = large_value1/DIST(vhelp,posA+nrmonA-1);
    mon2->x = vhelp->x + large_value2 * ((posA+nrmonA-1)->x - vhelp->x);
    mon2->y = vhelp->y + large_value2 * ((posA+nrmonA-1)->y - vhelp->y);
    mon2->z = vhelp->z + large_value2 * ((posA+nrmonA-1)->z - vhelp->z);
 
 
    // copy first
    VCOPY(mon1,posA);
    (posA)->monomer_number=1;
 
    // copy chain
    j=1;
    for(i=0;i<nrmonA;i++) {
      VCOPY((knotB+i),posA+j);
      ((posA+j)->monomer_number)=j+1;
      j++;
    }
 
    // copy last
    VCOPY(mon2,posA+nrmonA+1);
    (posA+nrmonA+1)->monomer_number=nrmonA+2;
 
    nrmonA+=2;
 
    for(i=0;i<nrmonA;i++)
     VPRINT(posA+i);
}

