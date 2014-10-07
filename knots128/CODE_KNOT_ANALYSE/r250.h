/***************************************************/
/*                                                 */
/* File: r250.h                                    */
/*                                                 */
/* Prototypes for the R250 random number generator */
/*                                                 */
/***************************************************/


/*#define PROBLEM_R250*/
/*#define SERIOUS_R250*/

#ifndef R250_H
#define R250_H

static char
r250_h[]="$Id: r250.h,v 1.5 1991/05/15 11:23:19 " 
         "reger Exp reger $";

double   double_r250( void );
float    float_r250( void );
void     double_r250_vector( register double *x, int n );
void     float_r250_vector( register float *x, int n );
float    float_rand( void );
double   gauss_r250( double mean, double variance );
int      int_rand( void );
int      r250_random( void );
void     r250_restart( int table[] );
void     r250_save( int table[] );
void     r250_srandom( unsigned int seed );
void     r250_vector( register int *x, int n );
int      rand_save( void );
/*void   srand( int seed );*/

#endif
