/***********************************************************
 * alexander_polynomial.c: determines alexander polynomial *
 *                                                         *
 * last modification: 09/21/2004                           *
 ***********************************************************/


#include "element.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int     mcs, polA, nrmonA;
extern MYVEC   *posA;
extern int     nrmon_knotA, nrmon_knotB;
extern MYVEC   *knotA, *knotB;


/******************************************************************
 * alexander_polynomial(double t)                                 *
 *           determines alexander polynomial at specified value t *
 *                                                                *
 * double t: value at which alexander polynomial is calculated    *
 *           coordinates of loop are taken from global matrix     *
 *           knotB                                                *
 *                                                                *
 * returns:  value of alexander polynomial at t                   *
 *                                                                *
 * last modification: 09/21/04                                    *
 ******************************************************************/
                                                                                                                             
double alexander_polynomial (double t)
{

  int    i,j,list_ctr=0,region_ctr=0;
  double a1, a2, m1, m2, x, y, z1, z2, z3, r, s;
  double distance=0;
  double dist[10*NMONOMAXA];
  INTERSECTION_LIST_ELEMENT *list, *intersection_list;

  double d,determinant=0;
  long double **alexander_matrix;
  int    indx[10*NMONOMAXA];

  /* allocate memory for list of intersecting bonds */
  list = (INTERSECTION_LIST_ELEMENT *) calloc(10*NMONOMAXA, sizeof(INTERSECTION_LIST_ELEMENT));
  intersection_list = list;

  /*----------------------------------------------------------------------------*
   | 1. project n bonds into xy-plane - create list of N intersecting bonds:    |
   |                                                                            |
   |    output: list of N crossings with                                        |
   |            intersection_list->section_i,j,k                                |
   |               regions assigned to the overpass (i) and the underpass (j,k) |
   |               of each crossing                                             |  
   |            intersection_list->right_handed (YES or NO)                     |
   |               determines whether or not crossing is right-handed           |   
   |                                                                            |
   |    (currently O(n^2) - could be improved to O(n) )                         |
   *----------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------------------------*
   | 1.1. calculate distances (along the polymer) between first monomer and successive monomers |
   |      (the resulting array will be used in section 1.2.2.1.)                                |
   *--------------------------------------------------------------------------------------------*/

   /*  printf("!%d\n",knotB->monomer_number);*/ 
   dist[(knotB)->monomer_number]=0; 
   // dist[1]=0; 
   for(i=0;i<nrmon_knotB-1;i++) {
     distance += sqrt( ((knotB+i)->x - (knotB+i+1)->x) * ((knotB+i)->x - (knotB+i+1)->x)  
                     + ((knotB+i)->y - (knotB+i+1)->y) * ((knotB+i)->y - (knotB+i+1)->y) );
     dist[(knotB+i+1)->monomer_number]=distance;  /* monomer_number 2,3,4, ..., n */   
   }

  /*----------------------------------------------------------------------*
   | 1.2. loop over all bonds (potential crossings) in 2D (xy) projection |
   *----------------------------------------------------------------------*/

   for(i=0;i<nrmon_knotB;i++) { 
     for(j=i+2;j<nrmon_knotB-1;j++) {   

     /* printf("%d %d\n",(knotB+i)->monomer_number,(knotB+j)->monomer_number); */ 
     
     /*-------------------------------------------------------------------------------------------*
      | 1.2.1. determine intersection point of lines through two specified bonds in xy-projection |   
      *-------------------------------------------------------------------------------------------*/
 
       m1 = ( (knotB+i)->y - (knotB+i+1)->y ) / ( (knotB+i)->x - (knotB+i+1)->x ); /* y1 = a1 + m1*x1 */
       m2 = ( (knotB+j)->y - (knotB+j+1)->y ) / ( (knotB+j)->x - (knotB+j+1)->x ); /* y2 = a2 + m2*x2 */
       
       a1 = (knotB+i)->y - m1 * (knotB+i)->x;
       a2 = (knotB+j)->y - m2 * (knotB+j)->x;

       x = (a1-a2) / (m2-m1);  /* x and y are coordinates of intersection point */
       y = a1 + m1 * x; 

     /*-----------------------------------------------------*
      | 1.2.2. check whether or not the two bonds intersect |
      *-----------------------------------------------------*/

       r = (x - (knotB+i)->x) / (-(knotB+i)->x + (knotB+i+1)->x);  
       s = (x - (knotB+j)->x) / (-(knotB+j)->x + (knotB+j+1)->x); 

       if( r>=0 && r<=1 && s>=0 && s<=1) {  /* YES? -> */
          
         /* determine z-coordinates of 2d intersection point on both (3d) bonds -> overpass/underpass */
         z1 = (knotB+i)->z + r * (-(knotB+i)->z + (knotB+i+1)->z);
         z2 = (knotB+j)->z + s * (-(knotB+j)->z + (knotB+j+1)->z);

        /*----------------------------------------------------------------------------*
         | 1.2.2.1. calculate distance (along the polymer) between first monomer and  |
         |          underpass (distance1) and first monomer and overpass (distance2), |
         |          assign values to interaction_list                                 |
         *----------------------------------------------------------------------------*/

         if(z1>z2) { /* first bond is overpassing second */
 
         /* printf("overpass\n"); */
        
         /* calculate distance of underpass to first monomer */
         distance = sqrt( ((knotB+j)->x - x) * ((knotB+j)->x - x) + ((knotB+j)->y - y) * ((knotB+j)->y - y) )
                    + dist[(knotB+j)->monomer_number];
         intersection_list->distance1 = distance;   

         /* calculate distance of overpass to first monomer */
         distance = sqrt( ((knotB+i)->x - x) * ((knotB+i)->x - x) + ((knotB+i)->y - y) * ((knotB+i)->y - y) )
                    + dist[(knotB+i)->monomer_number];
         intersection_list->distance2 = distance;
         }
       
         else {      /* first bond is underpassing second */

	 /* printf("underpass\n"); */

         /* calculate distance of underpass to first monomer */
         distance = sqrt( ((knotB+i)->x - x) * ((knotB+i)->x - x) + ((knotB+i)->y - y) * ((knotB+i)->y - y) )
                    + dist[(knotB+i)->monomer_number];
         intersection_list->distance1 = distance;
                                                               
         /* calculate distance of overpass to first monomer */                                                              
         distance = sqrt( ((knotB+j)->x - x) * ((knotB+j)->x - x) + ((knotB+j)->y - y) * ((knotB+j)->y - y) )
                    + dist[(knotB+j)->monomer_number];
         intersection_list->distance2 = distance;
         }

         /* printf("%lf %lf\n",intersection_list->distance1,intersection_list->distance2); */ 

        /*-------------------------------------------------------------*
         | 1.2.2.2. determine whether or not crossing is right-handed, |
         |          assign intersection_list->right_handed             |
         *-------------------------------------------------------------*/

         if(z1>z2) { /* first bond is overpassing 2nd */

           /* z3 = z-component of cross product between both intersecting bonds */
           /* right-handed: z3>0, left-handed: z3<0                             */
           z3 = (-(knotB+i)->x + (knotB+i+1)->x) * (-(knotB+j)->y + (knotB+j+1)->y) 
              - (-(knotB+i)->y + (knotB+i+1)->y) * (-(knotB+j)->x + (knotB+j+1)->x);
         }

         else {      /* first bond is underpassing 2nd */
           z3 = (-(knotB+j)->x + (knotB+j+1)->x) * (-(knotB+i)->y + (knotB+i+1)->y)
              - (-(knotB+j)->y + (knotB+j+1)->y) * (-(knotB+i)->x + (knotB+i+1)->x);
         }

         if(z3>0)
           intersection_list->right_handed = YES;
         else
           intersection_list->right_handed = NO;

         /* printf("!%d %d\n",(knotB+i)->monomer_number,(knotB+j)->monomer_number); */  

         intersection_list ++;
         list_ctr++;

       } /* end if(r>=0 && r<=1 && s>=0 && s<=1)  */
     } /* end for(j) */ 
   } /* end for(i) */

  if(list_ctr<3) {
    free(list);
    return(1);
  }

  /*---------------------------------------------------*
   | 1.3. assign regions to intersection list:         |
   |      region l is defined from underpasses (l-1,l] |
   |                                                   |
   |      (currently O(N^2), N = number of crossings   |
   |       - could be improved to O(NlogN))            |
   *---------------------------------------------------*/

   /* determine regions of underpassing bonds j,k */
   intersection_list=list;
   for(i=0;i<list_ctr;i++) {
     region_ctr=0;
     for(j=0;j<list_ctr;j++) {
       if( (i!=j) && ((list+i)->distance1 > (list+j)->distance1) ) {
         region_ctr++;
         /* printf("%d %d %lf %lf\n",i,j,(list+i)->distance1,(list+j)->distance1); */
       }
     }
     intersection_list->section_j=((region_ctr)%list_ctr)+1;
     intersection_list->section_k=((region_ctr+1)%list_ctr)+1;  /* k=j+1 ! */
     intersection_list++;
   } 

   /* determine regions of overpassing bond i */
   intersection_list=list;
   for(i=0;i<list_ctr;i++) {
     region_ctr=0;
     for(j=0;j<list_ctr;j++) {
       if( ((list+i)->distance2 > (list+j)->distance1) ) {
         region_ctr++;
       }      
     }
     intersection_list->section_i=((region_ctr)%list_ctr+1); 
     intersection_list++;
   }

   /* for(i=0;i<list_ctr;i++) {
   printf("%d %d %d %d\n",(list+i)->right_handed,(list+i)->section_i,(list+i)->section_j,(list+i)->section_k);
   } */


  /*----------------------------------------*
   | 2. assign elements to alexander matrix |
   *----------------------------------------*/

   /* allocate memory for alexander matrix (from Numerical Recipies) */
   alexander_matrix = dmatrix(1,list_ctr-1,1,list_ctr-1);
                        
   /* initialize matrix */                                            
   for(i=1;i<=list_ctr-1;i++)
     for(j=1;j<=list_ctr-1;j++) {
        alexander_matrix[i][j]=0;
   }

   /* assign elements - last row and last column are omitted */
   intersection_list=list;
   for(i=1;i<=list_ctr-1;i++) {
     if (intersection_list->right_handed == YES) {                               /* crossing right-handed */

       if (intersection_list->section_i  == intersection_list->section_j) {      /* i=j */
          if (intersection_list->section_i != list_ctr)
             alexander_matrix[i][intersection_list->section_i]=-t;
          if (intersection_list->section_k != list_ctr)
             alexander_matrix[i][intersection_list->section_k]=t;
       }

       else if (intersection_list->section_i  == intersection_list->section_k) { /* i=k=j+1 */
          if (intersection_list->section_i != list_ctr)
             alexander_matrix[i][intersection_list->section_i]=1;
          if (intersection_list->section_j != list_ctr)
             alexander_matrix[i][intersection_list->section_j]=-1;
       }

       else  {
          if (intersection_list->section_i != list_ctr)                          /* i!=j && i!=j+1 */
             alexander_matrix[i][intersection_list->section_i]=1-t;
          if (intersection_list->section_j != list_ctr)
             alexander_matrix[i][intersection_list->section_j]=-1;
          if (intersection_list->section_k != list_ctr)
             alexander_matrix[i][intersection_list->section_k]=t;
       }
     } /* end if */
  

     if (intersection_list->right_handed == NO) {                                /* crossing left-handed */
       if (intersection_list->section_i  == intersection_list->section_j) {      /* i=j */
          if (intersection_list->section_i != list_ctr)
             alexander_matrix[i][intersection_list->section_i]=1;
          if (intersection_list->section_k != list_ctr)
             alexander_matrix[i][intersection_list->section_k]=-1;
       }
                                                                                                                             
       else if (intersection_list->section_i  == intersection_list->section_k) { /* i=k=j+1 */
          if (intersection_list->section_i != list_ctr)
             alexander_matrix[i][intersection_list->section_i]=-t;
          if (intersection_list->section_j != list_ctr)
             alexander_matrix[i][intersection_list->section_j]=t;
       }
                                                                                                                             
       else  {
          if (intersection_list->section_i != list_ctr)                          /* i!=j && i!=j+1 */
             alexander_matrix[i][intersection_list->section_i]=1-t;
          if (intersection_list->section_j != list_ctr)
             alexander_matrix[i][intersection_list->section_j]=t;
          if (intersection_list->section_k != list_ctr)
             alexander_matrix[i][intersection_list->section_k]=-1;
       }
     } /* end if */
     intersection_list++;
   } /* end for */

  /*--------------------------------------------------*
   | 3. calculate determinant (=alexander polynomial) | 
   |    (LU decomposition -> O(N^3))                  |  
   |                                                  |
   |    from Numerical Recipies                       |
   *--------------------------------------------------*/

  ludcmp(alexander_matrix,list_ctr-1,indx,&d); 

  for(i=1;i<=list_ctr-1;i++) {
     d*=alexander_matrix[i][i];
  }
 
  determinant = d;

  free_dmatrix(alexander_matrix,1,list_ctr-1,1,list_ctr-1); 

  free(list);
  return(sqrt(determinant*determinant));
}
