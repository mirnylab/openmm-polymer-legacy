/*******************************************************************
 * determine_knot_size.c: functions to determine the size of knots *
 *                                                                 *
 * last modification: 02/15/2005                                   *
 *******************************************************************/


#include "element.h"
#include "r250.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int     mcs, polA, nrmonA;
extern MYVEC   *posA;
extern int     nrmon_knotA, nrmon_knotB;
extern MYVEC   *knotA, *knotB;
extern DATA      data;

extern int     distance_to_end;

extern char    *redknotout;
extern int     redknotcnt;


/*****************************************************************
 * test_for_knots: - copies relevant section of posA to knotB    *
 *                 - connects ends of open chain before or       *
 *                   after Taylor reduction                      *
 *                 - calculates Alexander polynomial             *
 *                                                               *
 * input: mode 0  : connects ends before Taylor reduction        *
 *             1  :    "     "    after     "       "            *
 *             2  : Taylor reduction without closure (loop)      *
 *   closure_mode : see function definition                      *
 *                  (only 0 and 2 are implemented so far)        *
 *             t  : position at which Alexander polynomial is    *
 *                  evaluated                                    *
 *                                                               *
 * The first and the last particle are always part of the chain. *
 * In addition, the section between "start" and "end" (monomer   *
 * positions in posA) is included. r_start (r_end) demark the    *
 * first (last) monomer after (before) the first (last) particle *
 * in the chain after the Taylor reduction.                      *
 *                                                               *
 * returns: |Alexander polynomial(t)*Alexander polynomial(1/t)|  * 
 *                                                               *
 * last modification: 02/15/2005                                 *
 *****************************************************************/


double test_for_knots (int mode, int closure_mode, double t, int start, int end, int *r_start, int *r_end)
{
  int i;

  if(mode==0) {
     copy_MYVEC_array_3(posA,knotB,start,end);
     nrmon_knotB=2+end-start+1;
     closure(closure_mode);

     // for(i=0;i<nrmon_knotB;i++)
     //  VPRINT(knotB+i);

     reduce_complexity_2();
     *r_start=((knotB+1)->monomer_number)-1;  /* positions in posA and monomer_number differ by 1 */
     *r_end  =((knotB+nrmon_knotB-2)->monomer_number)-1;
  }

  if(mode==1) {
     copy_MYVEC_array_3(posA,knotB,start,end);
     nrmon_knotB=2+end-start+1;

     /* for(i=0;i<nrmon_knotB;i++)
     VPRINT(knotB+i);
     printf("\n"); */

     reduce_complexity_2();
     *r_start=((knotB+1)->monomer_number)-1;
     *r_end  =((knotB+nrmon_knotB-2)->monomer_number)-1;
     closure(closure_mode); 
  }

  if(mode==2) {
     copy_MYVEC_array_3(posA,knotB,start,end);
     nrmon_knotB=2+end-start+1;
     reduce_complexity_2();
     *r_start=((knotB+1)->monomer_number)-1;
     *r_end  =((knotB+nrmon_knotB-2)->monomer_number)-1;
  }

  return(alexander_polynomial(t)*alexander_polynomial((1/t)));  
}


/**************************************************
 * generate_educated_guess_for_knot_size: dito    *
 *                                                *
 * input: p: |A(t)*A(1-t)| of knot type we are    *
 *        looking for                             *                      
 *        start: initial guess for first monomer  *
 *               after first (compare with above) *
 *        end: " " " last " before last           *
 *                                                *
 * r_start: e.g. for first monomer after first    *
 * r_end  : " "   "  last     "    before last    *
 *                                                *
 * last modification: 02/15/2005                  *
 **************************************************/
   
int generate_educated_guess_for_knot_size(double p, int start, int end, int *r_start, int *r_end) {

  int    s=start, e=end, ctr=0;
  int    r_s=start, r_e=end;
  int    s_old, e_old;
  int    temp1, temp2;
  double polynomial=p,polynomial_old;
 
  /*---------------------------------------------*
   | apply Taylor reduction until unknot appears |
   | 1. from the beginning of the chain          |     
   *---------------------------------------------*/   

  while(polynomial-1>0.000001) {   /* search until unknot appears */
                   
    s_old = s;
    s     = r_s;  
    
    polynomial=test_for_knots (1, 0, -1.1, s, e, &r_s, &r_e);
    ctr++;
    /* printf("0 %lf %d %d\n",polynomial,s,e); */ 

    while(r_s==s && (polynomial-1>0.000001) ) {   /* stuck */
      s_old=s;
      s++;
      polynomial=test_for_knots (1, 0, -1.1, s, e, &r_s, &r_e);      
      ctr++;
      /* printf("1 %lf %d %d\n",polynomial,s,e); */
    }
  }

 /*---------------------------------------------------*
  | add until resulting polymer section contains knot |
  *---------------------------------------------------*/
                                                                                                                             
  while(sqrt((polynomial-p)*(polynomial-p)) > 0.000001) {  /* while(polynomial != p) */
    if(s>1)  s=s-1;
    else     s=1;
                                                                                                                             
    polynomial=test_for_knots (1, 0, -1.1, s, e, &r_s, &r_e);
    ctr++;
    /* printf("2 %lf %d %d\n",polynomial,s,e); */
  }

  /*---------------------------------------------*
   | 2. from the end of the chain                |
   *---------------------------------------------*/

  r_e=end;

  while(polynomial-1>0.000001) {   /* search until unknot appears */
                                                                                                                             
    e     = r_e;
    polynomial=test_for_knots (1, 0, -1.1, s, e, &r_s, &r_e);
    ctr++;            
                                                                                                              
    while(r_e==e && (polynomial-1>0.000001) ) {   /* stuck */
      e--;
      polynomial=test_for_knots (1, 0, -1.1, s, e, &r_s, &r_e);
      ctr++;
    }
  } 

  /* printf("%d\t%d\t%d\t",ctr,s,e); */

  /*---------------------------------------------------*
   | add until resulting polymer section contains knot |
   *---------------------------------------------------*/

  while(sqrt((polynomial-p)*(polynomial-p)) > 0.000001) {  /* while(polynomial != p) */

    if(e<nrmonA-2-1)  e=e+1;
    else              e=nrmonA-2;

    polynomial=test_for_knots (1, 0, -1.1, s, e, &r_s, &r_e);
  }

  *r_start=s;
  *r_end  =e; 

  // printf("%d\t%d\t",s,e); 
  return(OK);
}


/***********************************************************************************
 * determine_knot_size: dito                                                       *
 *                                                                                 *
 * input:   p value for the product of the alexander polynomials (A(t1)*A(1/t1))   *
 *          s defines first monomer to be included in the middle part of the chain *
 *          e defines last monomer to be included in the middle part of the chain  *
 * returns: knot size (in monomer units)                                           *
 *                                                                                 *
 * last modification: 02/15/2005                                                   *
 ***********************************************************************************/
  
int determine_knot_size(double p, int s, int e, int flag)
{
  double polynomial;
  int i=0;
  int r_e,r_s;
  double t1,t2;
  int start=0,end=0,size;
  //MYVEC mon[1];
    
  char filenameout[30]="conf_red";
  
  /*------------------------------------------------------------------------*
   | The first and the last particle are always part of the chain.          |
   | Apart from that the rest of the chain consists of the section          |
   | between (and including) s and e (monomer positions in posA).           |
   | The programs starts with one end of the remaining chain and            |
   | removes the monomer next to the end particle. Then it applies Taylor's |
   | reduction and calculates the product of the Alexander                  |
   | polynomials at value t1 and 1/t1 (A(t1)*A(1/t1)).                      |
   | If the structure stays the same, variable "start" (monomer position in |
   | posA) is set to the monomer number which has been removed.             |
   | The iteration stops when the structure becomes an unknot.              |
   | The proceedure is repeated from the other side of the chain.           |
   | This time variable "end" is determined.                                |
   | The size of the knot is defined as "end"-"start" (number of bonds!)    |
   *------------------------------------------------------------------------*/
  
  t1=-1.1;      /* t1 and t2 should be the same as in check_for_knots */
  t2=1/t1;
  polynomial=p;

  start=s;
  end=e;

 
  /*--------------------------------*
   | 1. reduce chain from beginning |
   *--------------------------------*/
  
  for(i=s;polynomial-1>0.000001;i++){      /* stop when unknot is reached */

    polynomial=test_for_knots (1, 0, -1.1, i, e, &r_s, &r_e);
    /* printf("%lf\n", polynomial); */
    if( sqrt((polynomial-p)*(polynomial-p))  < 0.000001) {   /* set variable "start" if "polynomial==p" */
      start=i;
    }
  }
  polynomial=p;
  
  /*--------------------------*
   | 2. reduce chain from end |
   *--------------------------*/
  
  for(i=e;polynomial-1>0.000001;i--){

    polynomial=test_for_knots (1, 0, -1.1, start, i, &r_s, &r_e);
    /* printf("%lf\n", polynomial); */
    if( sqrt((polynomial-p)*(polynomial-p)) < 0.000001) {
      end=i;
    }
  }
 
  size=end-start; /* number of bonds ! */

  if(flag==0)  
    printf("%d\t%d\t",start,end);
    //printf("%d\t\t\t\t%d\t",start,end);

  polynomial=test_for_knots (1, 0, -1.1, start, end, &r_s, &r_e);

  /* prepare conf_red */
  /* knotB->x/=1000000;
  knotB->y/=1000000;
  knotB->z/=1000000;
  (knotB+nrmon_knotB-2)->x/=1000000;
  (knotB+nrmon_knotB-2)->y/=1000000;
  (knotB+nrmon_knotB-2)->z/=1000000;  */

  knotB->x= (knotB+1)->x + 0.000001 * ( -((knotB+1)->x) + (knotB->x) );
  knotB->y= (knotB+1)->y + 0.000001 * ( -((knotB+1)->y) + (knotB->y) );
  knotB->z= (knotB+1)->z + 0.000001 * ( -((knotB+1)->z) + (knotB->z) );

  (knotB+nrmon_knotB-2)->x= (knotB+nrmon_knotB-3)->x + 0.000001 * ( -(knotB+nrmon_knotB-3)->x + ((knotB+nrmon_knotB-2)->x) );
  (knotB+nrmon_knotB-2)->y= (knotB+nrmon_knotB-3)->y + 0.000001 * ( -(knotB+nrmon_knotB-3)->y + ((knotB+nrmon_knotB-2)->y) );
  (knotB+nrmon_knotB-2)->z= (knotB+nrmon_knotB-3)->z + 0.000001 * ( -(knotB+nrmon_knotB-3)->z + ((knotB+nrmon_knotB-2)->z) );



  if(redknotout){
//      char *p = NULL;
//      asprintf(&p,"%s.%d", redknotout, redknotcnt++);
//      if(p){
//	  sysout(p,knotB,nrmon_knotB-1,mcs);
//	  free(p);
//      }

      if(!redknotcnt++) // XXX
	  sysout(redknotout,knotB,nrmon_knotB-1,mcs);
  }




  return(size);
}


/***************************************************************************
 * prepare_loop_for_knot_analysis(): check if monomer 0 is part of largest *
 *                                   and second largest bond to ensure     *
 *                                   that it is not part of the knot       *
 *                                                                         *
 * last modification: 04/05/2005                                           *
 ***************************************************************************/

int prepare_loop_for_knot_analysis()
{
  int i,j,ctr=0;
  int start,end;
  int bondlength=0,distance,longest_bond;
  
  /*---------------------------------------------------* 
   | 1. Check if monomer 0 is part of the longest bond |       
   *---------------------------------------------------*/

  for(i=0;i<nrmon_knotB-1;i++){
  
      distance = (knotB+i+1)->monomer_number - (knotB+i)->monomer_number; 

      if(distance > bondlength){
        start      = ((knotB+i)   -> monomer_number)-1; /* start and end refer to positions in posA array */
        end        = ((knotB+i+1) -> monomer_number)-1;
        bondlength = distance;
      }
  }

  /* printf("%d\t%d\t",start,end); */

  if(start==0 || end==nrmonA-1)
    {}                         /* Monomer 0 is part of the longest bond ! -> continue */
  else
    return(NO);                /* No? -> stop function */


  /*---------------------------------------------------------------*
   | 2. Check if monomer 0 is also part of the second longest bond |
   |    (nearly identical to above)                                |
   *---------------------------------------------------------------*/

  longest_bond=bondlength;
  bondlength=0;
  for(i=0;i<nrmon_knotB-1;i++){
                                                                                                                            
      distance = (knotB+i+1)->monomer_number - (knotB+i)->monomer_number;
                                                                                                                             
      if( (distance > bondlength)  && (distance!=longest_bond) ){
        start      = ((knotB+i)   -> monomer_number)-1; 
        end        = ((knotB+i+1) -> monomer_number)-1;
        bondlength = distance;
      }
  }
  
  /* printf("%d\t%d\t",start,end);*/

  if(start==0 || end==nrmonA-1)
    return(YES);                       /* Monomer 0 is also part of the second longest bond ! */
  else
    return(NO);
}


/*********************************************************
 * change_position_of_first_monomer()                    *
 *    draw a random number from 0 to nrmonA-2            *
 *    shift posA such that drawn monomer_number becomes  *
 *    first monomer in loop                              *
 *                                                       *  
 * last modification: 04/05/2005                         *
 *********************************************************/

void change_position_of_first_monomer(int position){

int new_first_position;
int i,j;

  /*-----------------------*
   | 1. draw random number |
   *-----------------------*/

  new_first_position=position;

  /* new_first_position=(int) ( double_r250() * (nrmonA-2) ) +1;
  if(new_first_position==(nrmonA-1)) {new_first_position--;} */

  /* printf("%d\n",new_first_position); */

  /*---------------*
   | 2. shift loop |
   *---------------*/

  for(i=new_first_position,j=0;i<nrmonA-1;i++) {      
    VCOPY(posA+NMONOMAXA+10+i,posA+j);                           
    (posA+j)->monomer_number=j+1;
    /* printf("%d %d\n",i,(posA+j)->monomer_number); */
    j++;
  }
  for(i=0;i<=new_first_position;i++) {
    VCOPY(posA+NMONOMAXA+10+i,posA+j);
    (posA+j)->monomer_number=j+1;
    /* printf("%d %d\n",i,(posA+j)->monomer_number);*/
    j++;
  }
  nrmonA=j;
  /* printf("\n%d\n\n",nrmonA); */
  
}



/***************************************************************
 * copy_MYVEC_array: copies array_a of struct MYVEC to array_b *
 *                   starting at position "start" and ending   *
 *                   at position "end"                         *
 *                                                             *
 * last modification: 09/29/2004                               *
 ***************************************************************/

void copy_MYVEC_array (MYVEC *array_A, MYVEC *array_B, int start, int end)
{
int i,j=0;

  for(i=start;i<end;i++,j++) {
    VCOPY(array_A+i,array_B+j);
    /* printf("%d\t",j);
    VPRINT(array_B+j); */
  }
}

void copy_MYVEC_array_2 (MYVEC *array_A, MYVEC *array_B, int start, int end)
{
int i,j=0;
                       
  VCOPY(posA,array_B);                /* copy first monomer from posA to first position in array_B */
  j++;  
                                                                                                    
  for(i=start;i<end;i++,j++) {        /* copy designated monomers from the middle */
    VCOPY(array_A+i,array_B+j);       /* start,i demark position in array - not monomer_number!  */
  }                                   /* last position(!) to be copied: end-1 */

  VCOPY(posA+nrmonA-1,array_B+j);     /* copy last monomer from posA to las position in array_B */ 
  (array_B+j)->monomer_number=(posA+nrmonA-1)->monomer_number; 
}

void copy_MYVEC_array_3 (MYVEC *array_A, MYVEC *array_B, int start, int end)
{
int i,j=0;
                                                                                                                             
  VCOPY(posA,array_B);                /* copy first monomer from posA to first position in array_B */
  j++;
                                                                                                                             
  for(i=start;i<=end;i++,j++) {       /* copy designated monomers from the middle */
    VCOPY(array_A+i,array_B+j);       /* start,end,i,j demark position in array - not monomer_number!  */
  }                                  
                                                                                                                             
  VCOPY(posA+nrmonA-1,array_B+j);     /* copy last monomer from posA to last position in array_B */
  (array_B+j)->monomer_number=(posA+nrmonA-1)->monomer_number;
}
