/*************************************************************
 * element.h: header file                                    *
 *            contains substitution macros, function macros, *
 *            definitions of data structures and function    *
 *            prototypes                                     *
 *                                                           *
 * last modification : 04/19/2004                            *
 *************************************************************/


#define NMONOMAXA 50000

#define OK    0
#define ERROR 1
#define YES   0
#define NO    1


/*-----------------------------------*
 | macros for mathematical functions |
 *-----------------------------------*/

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#define ABS(x) (( (x)>0 )? (x) : -(x) )
#define MAX(x,y) (( (x)>(y) )? (x) : (y) )
#define MIN(x,y) (( (x)<(y) )? (x) : (y) )
#define SQUARE(x) ( (x)*(x) )
#define CUBE(x) ( (x)*(x)*(x) )
#define ASGN(x) (( (x)>0 )? (-1.0) : (1.0) )


/*------------------------------------------------* 
 | macros for mathematical functions with vectors |
 *------------------------------------------------*/ 

#define BETRAG(a)    (sqrt(SQUARE((a)->x)+SQUARE((a)->y)+\
                                          SQUARE((a)->z)))

#define BETRAG2(a)   (SQUARE((a)->x)+SQUARE((a)->y)+\
                                     SQUARE((a)->z))

#define DOT(a,b)     (((a)->x)*((b)->x) + ((a)->y) *\
((b)->y)+ ((a)->z) * ((b)->z))

#define CWINKEL(a,b) ((DOT(a,b))/(BETRAG(a)*BETRAG(b)))

#define DIST(a,b)    (sqrt(SQUARE((a)->x - (b)->x)+\
                           SQUARE((a)->y - (b)->y)+\
                           SQUARE((a)->z - (b)->z)))

#define DISTANCE(a,b) (sqrt(SQUARE((a)->x - (b)->x)+\
                            SQUARE((a)->y - (b)->y)+\
                            SQUARE((a)->z - (b)->z)))
#define DISTANCE2(a,b) (SQUARE((a)->x - (b)->x)+\
                        SQUARE((a)->y - (b)->y)+\
                        SQUARE((a)->z - (b)->z))

#define VCOPY(a,b)     (b)->monomer_number = (a)->monomer_number;\
                       (b)->x = (a)->x;\
                       (b)->y = (a)->y;\
                       (b)->z = (a)->z;

#define DIFF(a,b,c)    (c)->x = (a)->x - (b)->x;\
                       (c)->y = (a)->y - (b)->y;\
                       (c)->z = (a)->z - (b)->z

#define ADD(a,b,c)     (c)->x = (a)->x + (b)->x;\
                       (c)->y = (a)->y + (b)->y;\
                       (c)->z = (a)->z + (b)->z

#define SDOT(n,a,b)    (b)->x = (a)->x*(n);\
                       (b)->y = (a)->y*(n);\
                       (b)->z = (a)->z*(n)

#define VDOT(a,b,c)    (c)->x = (a)->y*\
                                (b)->z-(a)->z * (b)->y;\
                       (c)->y = (a)->z*\
                                (b)->x-(a)->x * (b)->z;\
                       (c)->z = (a)->x*\
                                (b)->y-(a)->y * (b)->x



#define VPRINT(a) printf("%d %f %f %f\n",(a)->monomer_number,(a)->x,(a)->y,(a)->z)
     

/*-----------------------*
 | structure definitions |
 *-----------------------*/

typedef struct o_element ELEMENT;
typedef struct o_box     BBOX;

typedef struct {
  int monomer_number;
  double x;
  double y;
  double z;
} MYVEC;

struct o_element {
  int monomer;               /* monomer number */
  struct o_element *before;  /* pointer to previous element  */
  struct o_element *after;   /* pointer to next element      */
};


/*                                             element 0     */
/*                                             after         */
/*                      element 1       <------before        */
/*                      after           ------>              */
/*  element 2     <-----before                               */
/*  after         ----->                                     */
/*  before                                                   */

struct o_box {
  struct o_element *first;  /* points somewhere in the list, */
                            /* there is no beginning         */
  int population;           /* number of list elements ,     */
                            /* number of particles in box    */
};

typedef struct DATA {
  double bondangleA;
  double bondlength1A;
  double bondlength2A;
  double bondlengthmaxA;
  double gyrA;
  double endA;
  double msdA;
} DATA;

typedef struct {
  double distance1;
  double distance2;
  int right_handed;
  int section_i;
  int section_j;
  int section_k;
} INTERSECTION_LIST_ELEMENT;


/*-----------------------------------------------------------*/

/*----------------------*
 | function declaration |
 *----------------------*/

/*------*
 | io.c |
 *------*/
     
/* reads monomer positions from conf_in */ 
int sysin(char *);         

/* write sim. specifications and monomer positions to "conf_out" */
int sysout(char *, MYVEC *, int, int);

/*---------------*
 | random_walk.c |
 *---------------*/

int generate_3d_random_walk();
int setup_loop();
int rotate_loop();

/*--------*
 | knot.c |
 *--------*/
                                                                                                                             
int check_for_knots();
void closure (int method);
void change_posA();

/*-----------------------*
 | determine_knot_size.c |
 *-----------------------*/

double test_for_knots (int mode, int closure_mode, double t, int start, int end, int *r_start, int *r_end); 
int generate_educated_guess_for_knot_size(double p, int start, int end, int *r_start, int *r_end);
int determine_knot_size(double p, int s, int e, int flag);
int prepare_loop_for_knot_analysis();
void change_position_of_first_monomer(int position);
void copy_MYVEC_array (MYVEC *array_A, MYVEC *array_B, int start, int end);
void copy_MYVEC_array_2 (MYVEC *array_A, MYVEC *array_B, int start, int end);
void copy_MYVEC_array_3 (MYVEC *array_A, MYVEC *array_B, int start, int end);

/*----------*
 | taylor.c |
 *----------*/
                                                                                                                             
/* reduces complexity of polymer by successively deleting redundant monomers */
int reduce_complexity();
int reduce_complexity_1();
int reduce_complexity_2();
                                                                                                                             
/* subfunction of reduce_complexity */
int check_crossing(MYVEC *A, MYVEC *B, MYVEC *C);
                                                                                                                             
/*------------------------*
 | alexander_polynomial.c |
 *------------------------*/
                                                                                                                             
double alexander_polynomial(double position);

/*----------------*
 | density_prof.c |
 *----------------*/

void calc_cm1(int pol, MYVEC *start, MYVEC *result);
void dens_prof(MYVEC *pos, int nrmon, double abstand);

/*------------*
 | analysis.c |
 *------------*/

void analysis(int zeit);
void bondwinkel();        /* calculates bondangles                      */
void bondlength();        /* bondlength, b**2, bmax                     */
void calc_polymer();      /* end-to-end distance and radius of gyration */
void calc_cm();           /* calculates center of mass                  */
void calc_msd();          /* calculates mean square displacement        */

/*-------------*
 | numerical.c |
 *-------------*/
  
/* functions from Numerical Recipies in C
   (used in alexander_polynomial.c to calculate the determinant of the alexander matrix) */
  
long double **dmatrix(int nrl,int nrh,int ncl,int nch);
void free_dmatrix(long double **m,int nrl,int nrh,int ncl,int nch);
void ludcmp(long double **a,int n,int *indx,double *d);
long double *vector(int nl,int nh);
void nrerror(char *error_text);
void free_vector(long double *v, int nl, int nh);
