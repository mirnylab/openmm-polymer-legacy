/***********************************************************
 * taylor.c: functions to delete redundant monomers from   *
 *           chain to reduce complexity before calculating *
 *           knot polynomials                              *
 *                                                         *
 * compare with Taylor WR, Nature 406, 916 (2000)          *       
 *                                                         *
 * last modification: 09/09/2004                           *
 ***********************************************************/


#include "element.h"
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

extern int     mcs, polA, nrmonA;
extern MYVEC   *posA;
extern int     nrmon_knotA, nrmon_knotB;
extern MYVEC   *knotA, *knotB;


/**************************************************************
 * reduce_complexity(): successively deletes elements in list *
 *                      if procedure doesn't destroy knots    *
 *                      (check_crossing()==NO, see below)     *
 *                                                            *
 * returns:		number of monomers left (2->no knot)  *
 * 			(>2 check with knot polynomials)      *
 *							      *		
 * last modification: 03/29/2004                              *
 **************************************************************/

int reduce_complexity ()
{
  int i, ctr=0;
  MYVEC *mon_knotA, *mon_knotB, *mon;
  //char filenameout[30]="conf_red";


 /*--------------------------*
  | copy posA-array to knotB |
  *--------------------------*/

  /* nrmon_knotB=nrmonA;
  mon=posA; 
  mon_knotB=knotB;

  for(i=0;i<nrmonA;i++){

    VCOPY(mon,mon_knotB);  
    
    mon++;
    mon_knotB++;    
  } */

  while(nrmon_knotB!=nrmon_knotA){ 

    /*-----------------------------------------------*
     | "copy" knotB to knotA by switching identities |
     *-----------------------------------------------*/

    mon=knotA;
    knotA=knotB;
    knotB=mon;
 
    nrmon_knotA=nrmon_knotB;
    nrmon_knotB=0;
  
    mon_knotA=knotA;
    mon_knotB=knotB;

    /*-----------------------------------------------------------------------------*
     | start at beginning of knotA and n=1                                         |
     |                                                                             |
     | |: check if removal of n+1th particle would interfere with any line segment |
     |      YES -copy nth and n+1th particle positions to knotB                    |
     |      NO  -copy nth particle to knotB, omit (delete) n+1th particle          |
     |    n++ :| (loop over all particles n in knotA)                              |
     *-----------------------------------------------------------------------------*/

    for(i=0; i<nrmon_knotA; i+=2){  /* i denotes how many steps mon_knotA advanced in knotA */
                                    /* f.i. i==nrmon_knotA-1 -> mon_knotA points on last    */
                                    /*      element in knotA */

      VCOPY(mon_knotA,mon_knotB);         /* copy 1st element of triangle */ 
      mon_knotA++;
      mon_knotB++; 
      nrmon_knotB++;
     
      if(i<nrmon_knotA-2){      /* at least 3 monomers must remain to form triangle ! */
        if(check_crossing(mon_knotA-1,mon_knotA,mon_knotA+1)==YES){  /* interference: */ 
          VCOPY((mon_knotA),(mon_knotB));         /* -> copy 2nd particle of triangle */  
	  mon_knotB++;                            /* else: omit (delete) 2nd particle */
	  nrmon_knotB++;
	}
	mon_knotA++;
      }
    } /* end for-loop */

    /*---------------------------------------------------------------* 
     | last particle of knotA also present in knot B ?               |
     |   YES: ok                                                     |
     |   NO:  copy last particle of knotA to last position of knot B |
     *---------------------------------------------------------------*/    

    if ( DISTANCE2(knotA+(nrmon_knotA-1),knotB+(nrmon_knotB-1)) != 0 ){ /* particles not identical ?  */
      VCOPY(knotA+(nrmon_knotA-1),knotB+(nrmon_knotB));                 /* be careful at this point ! */
      nrmon_knotB++;
    }     

    ctr++;    
    /* sprintf(filenameout,"k_%d",ctr);
    sysout(filenameout,knotB,nrmon_knotB,mcs); */ 

  } /* end while */

  /* sprintf(filenameout,"k_1");
  sysout(filenameout,knotB,nrmon_knotB,mcs); */ 

  /* printf("%d\t%d\t",(knotB+1)->monomer_number,(knotB+nrmon_knotB-2)->monomer_number); */
  return(nrmon_knotB); 
}


/**************************************************************
 * reduce_complexity_1(): as above but starts with end-monomer*
 *                                                            *
 * returns:		number of monomers left (2->no knot)  *
 * 			(>2 check with knot polynomials)      *
 *							      *		
 * last modification: 03/29/2004                              *
 **************************************************************/

int reduce_complexity_1 ()
{
  int i, ctr=0;
  MYVEC *mon_knotA, *mon_knotB, *mon;
  char filenameout[30]="conf_red";


 /*--------------------------*
  | copy posA-array to knotB |
  *--------------------------*/

  /* nrmon_knotB=nrmonA;
  mon=posA+(nrmonA-1); */ /* ! */
  /* mon_knotB=knotB;

  for(i=0;i<nrmonA;i++){

    VCOPY(mon,mon_knotB);  
    mon--;
    mon_knotB++;    
  } */

  while(nrmon_knotB!=nrmon_knotA){ 

    /*-----------------------------------------------*
     | "copy" knotB to knotA by switching identities |
     *-----------------------------------------------*/

    mon=knotA;
    knotA=knotB;
    knotB=mon;
 
    nrmon_knotA=nrmon_knotB;
    nrmon_knotB=0;
  
    mon_knotA=knotA;
    mon_knotB=knotB;

    /*-----------------------------------------------------------------------------*
     | start at beginning of knotA and n=1                                         |
     |                                                                             |
     | |: check if removal of n+1th particle would interfere with any line segment |
     |      YES -copy nth and n+1th particle positions to knotB                    |
     |      NO  -copy nth particle to knotB, omit (delete) n+1th particle          |
     |    n++ :| (loop over all particles n in knotA)                              |
     *-----------------------------------------------------------------------------*/

    for(i=0; i<nrmon_knotA; i+=2){  /* i denotes how many steps mon_knotA advanced in knotA */
                                    /* f.i. i==nrmon_knotA-1 -> mon_knotA points on last    */
                                    /*      element in knotA */

      VCOPY(mon_knotA,mon_knotB);         /* copy 1st element of triangle */ 
      mon_knotA++;
      mon_knotB++; 
      nrmon_knotB++;
     
      if(i<nrmon_knotA-2){      /* at least 3 monomers must remain to form triangle ! */
        if(check_crossing(mon_knotA-1,mon_knotA,mon_knotA+1)==YES){  /* interference: */ 
          VCOPY((mon_knotA),(mon_knotB));         /* -> copy 2nd particle of triangle */  
	  mon_knotB++;                            /* else: omit (delete) 2nd particle */
	  nrmon_knotB++;
	}
	mon_knotA++;
      }
    } /* end for-loop */

    /*---------------------------------------------------------------* 
     | last particle of knotA also present in knot B ?               |
     |   YES: ok                                                     |
     |   NO:  copy last particle of knotA to last position of knot B |
     *---------------------------------------------------------------*/    

    if ( DISTANCE2(knotA+(nrmon_knotA-1),knotB+(nrmon_knotB-1)) != 0 ){ /* particles not identical ?  */
      VCOPY(knotA+(nrmon_knotA-1),knotB+(nrmon_knotB));                 /* be careful at this point ! */
      nrmon_knotB++;
    }     

    ctr++;    
    /* sprintf(filenameout,"k_%d",ctr);
       sysout(filenameout,knotB,nrmon_knotB,mcs); */ 

  } /* end while */

  //sprintf(filenameout,"k_2");
  //sysout(filenameout,knotB,nrmon_knotB,mcs); 

  return(nrmon_knotB); 
}


/****************************************************************
 * reduce_complexity_2(): as above but alternates in directions *
 *                                                              *
 * returns:               number of monomers left (2->no knot)  *
 *                        (>2 check with knot polynomials)      *
 *                                                              *
 * last modification: 03/29/2004                                *
 ****************************************************************/
 
int reduce_complexity_2 ()
{
  int i, ctr=0, control_variable=0;
  MYVEC *mon_knotA, *mon_knotB, *mon;
  //char filenameout[30]="conf_red";
 
 
 /*--------------------------*
  | copy posA-array to knotB |
  *--------------------------*/
 
  /* nrmon_knotB=nrmonA;
  mon=posA;
  mon_knotB=knotB;
 
  for(i=0;i<nrmonA;i++){
 
    VCOPY(mon,mon_knotB); 
     
    mon++;
    mon_knotB++;
  } */
 
  while( (nrmon_knotB!=nrmon_knotA || control_variable==0) && nrmon_knotB!=2){
 
    if(nrmon_knotB==nrmon_knotA) {  /* this ensures that the chain is reduced from both directions */
      control_variable=1;           /* and doesn't get stuck if it cannot proceed in only one      */
    }
    else {
      control_variable=0;
    }

    /*-------------------------------------------------*
     | "copy" knotB to knotA - alternate in each round |
     *-------------------------------------------------*/
 
    nrmon_knotA=nrmon_knotB;
    mon=knotB+(nrmon_knotB-1); /* ! */
    mon_knotA=knotA;
  
    for(i=0;i<nrmon_knotB;i++){
  
      VCOPY(mon,mon_knotA);  /* copy coordinates from mon to mon_knotB */
      mon--;
      mon_knotA++;
    }
 
    nrmon_knotB=0;
   
    mon_knotA=knotA;
    mon_knotB=knotB;
 
    /*-----------------------------------------------------------------------------*
     | start at beginning of knotA and n=1                                         |
     |                                                                             |
     | |: check if removal of n+1th particle would interfere with any line segment |
     |      YES -copy nth and n+1th particle positions to knotB                    |
     |      NO  -copy nth particle to knotB, omit (delete) n+1th particle          |
     |    n++ :| (loop over all particles n in knotA)                              |
     *-----------------------------------------------------------------------------*/
 
    for(i=0; i<nrmon_knotA; i+=2){  /* i denotes how many steps mon_knotA advanced in knotA */
                                    /* f.i. i==nrmon_knotA-1 -> mon_knotA points on last    */
                                    /*      element in knotA */

      VCOPY(mon_knotA,mon_knotB);         /* copy 1st element of triangle */
      mon_knotA++;
      mon_knotB++;
      nrmon_knotB++;
 
      if(i<nrmon_knotA-2){      /* at least 3 monomers must remain to form triangle ! */
        if(check_crossing(mon_knotA-1,mon_knotA,mon_knotA+1)==YES){  /* interference: */
          VCOPY((mon_knotA),(mon_knotB));         /* -> copy 2nd particle of triangle */
          mon_knotB++;                            /* else: omit (delete) 2nd particle */
          nrmon_knotB++;
        }
        mon_knotA++;
      }
    } /* end for-loop */
 
    /*---------------------------------------------------------------*
     | last particle of knotA also present in knot B ?               |
     |   YES: ok                                                     |
     |   NO:  copy last particle of knotA to last position of knot B |
     *---------------------------------------------------------------*/

 
    if ( DISTANCE2(knotA+(nrmon_knotA-1),knotB+(nrmon_knotB-1)) != 0 ){ /* particles not identical ?  */
      VCOPY(knotA+(nrmon_knotA-1),knotB+(nrmon_knotB));                 /* be careful at this point ! */
      nrmon_knotB++;
    }
 
    /* printf("%d\n",nrmon_knotB);*/
    ctr++;
    /* sprintf(filenameout,"c_%d",ctr);
    sysout(filenameout,knotB,nrmon_knotB,mcs); */

  } /* end while */

  /* printf("%d\n",nrmon_knotB); */

  /* if order is switched (highest monomer_number first), switch back */
  /* if(knotB->monomer_number > (knotB+nrmon_knotB-1)->monomer_number) {
    nrmon_knotB=nrmon_knotA;
    mon=knotA+(nrmon_knotA-1); 
    mon_knotB=knotB;
                                                                                                                             
    for(i=0;i<nrmon_knotB;i++){
      VCOPY(mon,mon_knotB); 
      mon--;
      mon_knotB++;
    }
  } */
  
  /* if order is switched (highest monomer_number first), switch back */
  if(knotB->monomer_number > (knotB+nrmon_knotB-1)->monomer_number) {

    /* for(i=0;i<nrmon_knotB;i++)
      VPRINT(knotB+i); */

    for(i=0;i<nrmon_knotB;i++){
      VCOPY(knotB+nrmon_knotB-1-i,knotA+i);
    }
    nrmon_knotA=nrmon_knotB;
    
    for(i=0;i<nrmon_knotA;i++){
      VCOPY(knotA+i,knotB+i);
      /* VPRINT(knotB+i); */
    }  
  } 

  /* sprintf(filenameout,"c_1"); */
  /* if(nrmon_knotB!=2)
    printf("%d %d ",(knotB+1)->monomer_number,(knotB+nrmon_knotB-2)->monomer_number);
  else
    printf("0 0 "); */
 
  return(nrmon_knotB);
}


/*************************************************************
 * check_crossing(MYVEC *A,B,C): checks if triangle defined  *
 *                               by A, B and C is crossed by *
 *                               any line segment (i,i+1)    *
 *                               with i< A (in knotB) ||     *
 *                               i> C (in knotA)             *
 *                                                           *
 * remark:  A,B,C are pointers to vectors in knotA structure *
 *          vector B is checked to be removed !              *
 *                                                           *
 * returns: YES: crossing                                    *
 * 	    NO : no crossing                                 *
 *				                             *	
 * last modification: 04/08/2004                             *
 *************************************************************/

int check_crossing(MYVEC *A, MYVEC *B, MYVEC *C)
{
  double r,s,t;
  MYVEC *OA,*AB,*BC,*L1A,*L12,*mon_ctr;
  MYVEC *crossing_array;

  /*-----------------------------------------------------------------------*
   | 1.1 initialize array for storage of triangle and line segment vectors |
   *-----------------------------------------------------------------------*/

  crossing_array = (MYVEC *) calloc(5,sizeof(MYVEC)); 

  /*--------------------------------*
   | 1.2 assign vectors to triangle |
   *--------------------------------*/

  OA   = (crossing_array);
  AB   = (crossing_array+1);
  BC   = (crossing_array+2);
  L1A  = (crossing_array+3);
  L12  = (crossing_array+4);

  VCOPY( A, OA  );     /* A  = OA  */
  DIFF ( B, A , AB );  /* AB = B-A */
  DIFF ( C, B , BC );  /* BC = C-B */

  /*----------------------------------------------------------*
   | 2.1 loop over all line segments in knotB up to (A-2,A-1) |
   *----------------------------------------------------------*/

  for(mon_ctr=knotB; mon_ctr<knotB+(nrmon_knotB-2); mon_ctr++){  /* last monomer of knotB is A-1 ! */

    /*----------------------------*
     | 2.1.1 assign line segments |
     *----------------------------*/
    
    DIFF ( OA, mon_ctr, L1A);             /* L1A = OA-L1 */
    DIFF ( (mon_ctr+1), mon_ctr, L12);    /* L12 = L2-L1 */

    /*-------------------------------------------------------------------------------*
     | 2.1.2 determine intersection point of straight line through L12 and plane ABC |
     |       remark: the following section has been thorougly tested !               |
     *-------------------------------------------------------------------------------*/
   
    r =   ( - (L1A->y)*(L12->x)*(BC->z) - (L1A->x)*(L12->z)*(BC->y) 
            + (L1A->z)*(L12->x)*(BC->y) + (L1A->x)*(L12->y)*(BC->z) 
            + (L1A->y)*(L12->z)*(BC->x) - (L1A->z)*(L12->y)*(BC->x) )/ 
          ( - (L12->y)*(AB->x)*(BC->z)  + (L12->x)*(AB->y)*(BC->z) 
            - (L12->x)*(AB->z)*(BC->y)  - (L12->z)*(AB->y)*(BC->x) 
            + (L12->y)*(AB->z)*(BC->x)  + (L12->z)*(AB->x)*(BC->y) );

    s = - ( - (L1A->z)*(L12->y)*(AB->x) + (L1A->z)*(L12->x)*(AB->y) 
            + (L1A->y)*(L12->z)*(AB->x) - (L1A->x)*(L12->z)*(AB->y)
            - (L1A->y)*(L12->x)*(AB->z) + (L1A->x)*(L12->y)*(AB->z) ) /
          ( - (L12->y)*(AB->x)*(BC->z)  + (L12->x)*(AB->y)*(BC->z) 
            - (L12->x)*(AB->z)*(BC->y)  - (L12->z)*(AB->y)*(BC->x) 
            + (L12->y)*(AB->z)*(BC->x)  + (L12->z)*(AB->x)*(BC->y) );

    t =   ( - (L1A->y)*(AB->x)*(BC->z)  + (L1A->x)*(AB->y)*(BC->z) 
            + (L1A->y)*(AB->z)*(BC->x)  + (L1A->z)*(AB->x)*(BC->y)
            - (L1A->x)*(AB->z)*(BC->y)  - (L1A->z)*(AB->y)*(BC->x) ) / 
          ( - (L12->y)*(AB->x)*(BC->z)  + (L12->x)*(AB->y)*(BC->z) 
            - (L12->x)*(AB->z)*(BC->y)  - (L12->z)*(AB->y)*(BC->x)
            + (L12->y)*(AB->z)*(BC->x)  + (L12->z)*(AB->x)*(BC->y) ); 

    /*--------------------------------------------------------*
     | 2.1.3 check if intersection point lies within triangle |
     *--------------------------------------------------------*/
   
    if( (0<r) && (r<1) && (0<s) && (s<r) && (0<t) && (t<1) ) { /* modified: used to be <= ! */
      free(crossing_array);
      /* printf("%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",A->monomer_number,B->monomer_number,C->monomer_number,mon_ctr->monomer_number,(mon_ctr+1)->monomer_number,s,r,t); */
      return(YES);
      } 

  }  /* end for-loop (2.1) */

  /*--------------------------------------------------------------------------------------------*
   | 2.2 loop over all line segments in knotA starting with (C+1,C+2) up to last entry in knotA |
   |     apart from the for-loop this section is identical with 2.1                             |
   *--------------------------------------------------------------------------------------------*/

  for(mon_ctr=(C+1); mon_ctr<knotA+(nrmon_knotA-1); mon_ctr++){   /* first monomer is C+1 ! */

    /*----------------------------*
     | 2.2.1 assign line segments |
     *----------------------------*/

    DIFF ( OA, mon_ctr, L1A);             /* L1A = OA-L1 */
    DIFF ( (mon_ctr+1), mon_ctr, L12);    /* L12 = L2-L1 */

    /*-------------------------------------------------------------------------------*
     | 2.2.2 determine intersection point of straight line through L12 and plane ABC |
     *-------------------------------------------------------------------------------*/
    
    r =   ( - (L1A->y)*(L12->x)*(BC->z) - (L1A->x)*(L12->z)*(BC->y) 
            + (L1A->z)*(L12->x)*(BC->y) + (L1A->x)*(L12->y)*(BC->z) 
            + (L1A->y)*(L12->z)*(BC->x) - (L1A->z)*(L12->y)*(BC->x) )/ 
          ( - (L12->y)*(AB->x)*(BC->z)  + (L12->x)*(AB->y)*(BC->z) 
            - (L12->x)*(AB->z)*(BC->y)  - (L12->z)*(AB->y)*(BC->x) 
            + (L12->y)*(AB->z)*(BC->x)  + (L12->z)*(AB->x)*(BC->y) );

    s = - ( - (L1A->z)*(L12->y)*(AB->x) + (L1A->z)*(L12->x)*(AB->y) 
            + (L1A->y)*(L12->z)*(AB->x) - (L1A->x)*(L12->z)*(AB->y)
            - (L1A->y)*(L12->x)*(AB->z) + (L1A->x)*(L12->y)*(AB->z) ) /
          ( - (L12->y)*(AB->x)*(BC->z)  + (L12->x)*(AB->y)*(BC->z) 
            - (L12->x)*(AB->z)*(BC->y)  - (L12->z)*(AB->y)*(BC->x) 
            + (L12->y)*(AB->z)*(BC->x)  + (L12->z)*(AB->x)*(BC->y) );

    t =   ( - (L1A->y)*(AB->x)*(BC->z)  + (L1A->x)*(AB->y)*(BC->z) 
            + (L1A->y)*(AB->z)*(BC->x)  + (L1A->z)*(AB->x)*(BC->y)
            - (L1A->x)*(AB->z)*(BC->y)  - (L1A->z)*(AB->y)*(BC->x) ) / 
          ( - (L12->y)*(AB->x)*(BC->z)  + (L12->x)*(AB->y)*(BC->z) 
            - (L12->x)*(AB->z)*(BC->y)  - (L12->z)*(AB->y)*(BC->x)
            + (L12->y)*(AB->z)*(BC->x)  + (L12->z)*(AB->x)*(BC->y) ); 

    /*--------------------------------------------------------*
     | 2.2.3 check if intersection point lies within triangle |
     *--------------------------------------------------------*/

    if( (0<r) && (r<1) && (0<s) && (s<r) && (0<t) && (t<1) ) { /* modified: used to be <= ! */
      free(crossing_array);
      /* printf("%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",A->monomer_number,B->monomer_number,C->monomer_number,mon_ctr->monomer_number,(mon_ctr+1)->monomer_number,s,r,t); */
      return(YES);
      }  

  }  /* end for-loop (2.2) */

  free(crossing_array); 
  return(NO);
}
