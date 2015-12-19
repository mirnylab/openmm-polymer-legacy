/*********************************
 * io.c: input/output functions  *
 *                               *
 * last modification: 04/19/2004 *
 *********************************/


#include "element.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

extern int mcs, polA;
extern int nrmonA;
extern MYVEC *posA;
extern MYVEC *cm0A;
extern DATA  *data;

/********************************************************
 * sysin.c: reads monomer positions from file "conf_in" *
 *                                                      *
 * last modification: 04/19/2004                        *
 ********************************************************/

int sysin(char *name)
{
  int    i;
  MYVEC  *mon;
  FILE   *in;

  in = fopen(name,"r");

  /*-----------------------------------------------*
   | read simulation specifications from "conf_in" |
   *-----------------------------------------------*/ 

  /*-----------------*
   | allocate memory |
   *-----------------*/

  posA           = (MYVEC *) calloc(2*NMONOMAXA+20, sizeof(MYVEC));
  mon = posA;

  cm0A           = (MYVEC *) calloc(1, sizeof(MYVEC));
  data           = (DATA *) calloc(1,sizeof(DATA));
  
  if (posA==NULL) {
    fprintf(stderr,"ERROR: memory allocation denied"
                   "in function sysin\n");
    exit(ERROR);
  }

  /*-------------------------------------------*
   | read monomer positions from "conf_in" and |
   | update posA array                         |
   *-------------------------------------------*/

  fscanf(in,"t=%d\n\n", &mcs);
  fscanf(in,"%d\n", &polA);
  nrmonA = polA;
  for (i=0 ; i< polA ; i++){
    fscanf( in,"%d%lf%lf%lf",&(mon->monomer_number),&(mon->x),&(mon->y),&(mon->z));
    mon ++; 
  }

  polA=NMONOMAXA;

  fclose(in);
  return(0);

} /* end sysin */


/********************************************************
 * sysout.c: write monomer positions to file "conf_out" *
 *                                                      *
 * last modification: 04/08/2004                        *
 ********************************************************/

int sysout(char *name, MYVEC *pos_array, int nr_elements_in_pos_array, int mcs)
{
  int    i;
  MYVEC  *mon;                                                                             
  FILE   *out;
  
  if( (out = fopen(name,"w")) == NULL ){
    fprintf(stderr,"Can not open file %s for writing: %s\n", 
	name, strerror(errno));
    exit(1);
  }
  mon = pos_array;
 
  /*---------------------------------------*
   | write monomer positions to "conf_out" |
   *---------------------------------------*/

  fprintf(out,"t=%d\n\n",mcs); 
  /* fprintf(out,"%d\n\n", nr_elements_in_pos_array);*/
  fprintf(out,"%d\n", nr_elements_in_pos_array);
    for( i=0 ; i<nr_elements_in_pos_array; i++){
      /* fprintf(out,"%lf\t%lf\t%lf\n",mon->x, mon->y, mon->z); */
      fprintf(out,"%d\t%lf\t%lf\t%lf\n",mon->monomer_number,mon->x, mon->y, mon->z); 
      mon++ ; 
  }

  fclose(out);
  return(0);

  }  /* end sysout */
