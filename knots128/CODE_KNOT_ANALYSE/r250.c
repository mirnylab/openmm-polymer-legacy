/***************************************************/
/*                                                 */
/*  File: r250.c                                   */
/*                                                 */
/*  The R250 random number generator routines.     */
/*                                                 */
/***************************************************/

static char
rcs[] = "$Id: r250.c,v 1.1.1.1 2006/12/18 18:36:40 kolesov Exp $";

#include <math.h>
#include <stdio.h>
#include "r250.h"

#ifndef INT_MAX
#        define INT_MAX                  0x7FFFFFFF
#endif
#if INT_MAX > 0x7FFFFFFF
#        define INT_MAX                  0x7FFFFFFF
#endif

#define INVERSE_INT_MAX          (1.0/INT_MAX)

#define                TYPE_3                  3               
/* x**31 + x**3 + 1 */

#define                BREAK_3                128
#define                DEG_3                   31
#define                SEP_3                    3

#define P 250
#define Q 103

static int m[P]    =
  {  0x7be9c1bd, 0x088aa102, 0x3d38509b, 0x746b9fbe, 
     0x2d04417f, 0x775d4351, 0x53c48d96, 0x02b26e0b, 
     0x418fedcf, 0x19dbc19e, 0x78512adb, 0x1a1f5e2b, 
     0x307d7761, 0x6584c1f0, 0x24e3c36f, 0x2232310f, 
     0x2dac5ceb, 0x106e8b5a, 0x5a05a938, 0x5e6392b6,
     0x66b90348, 0x75264901, 0x4174f402, 0x618b18a4, 
     0x3a6ab4bf, 0x3c4bc289, 0x2657a9a9, 0x4e68589b, 
     0x09648aa6, 0x3fc489bb, 0x1c1b715c, 0x054e4c63, 
     0x484f2abd, 0x5953c1f8, 0x79b9ec22, 0x75536c3d, 
     0x50b10549, 0x4d7e79b8, 0x7805da48, 0x1240f318,
     0x675a3b56, 0x70570523, 0x2c605143, 0x17d7b2b7, 
     0x55dbc713, 0x514414b2, 0x3a09e3c7, 0x038823ff, 
     0x61b2a00c, 0x140f8cff, 0x61ebb6b5, 0x486ba354,
     0x0935d600, 0x2360aab8, 0x29f6bbf8, 0x43a08abf, 
     0x5fac6d41, 0x504e65a2, 0x1208e35b, 0x6910f7e7,
     0x1012ef5d, 0x2e2454b7, 0x6e5f444b, 0x58621a1b, 
     0x077816af, 0x6819306d, 0x4db58658, 0x58291bf8, 
     0x3597aa25, 0x45bb60a0, 0x6a6a0f11, 0x1cf1e57c, 
     0x361265c4, 0x16ca6054, 0x34c99833, 0x0bee2cd7, 
     0x680e7507, 0x6ed37bfa, 0x0f7650d6, 0x49c11513,
     0x02e308f9, 0x7162078c, 0x122cb868, 0x0c18defa, 
     0x14c2b244, 0x3c237460, 0x4fb969b9, 0x746f1f85, 
     0x0c71da02, 0x61c24d14, 0x5d80176d, 0x1c84c960, 
     0x0fe6a1cc, 0x4bdf5bb8, 0x74e6e37b, 0x175eb87b, 
     0x33f88c25, 0x429c69d3, 0x6f87d474, 0x6990364a,
     0x0857ca73, 0x59f1e385, 0x06821bc6, 0x3e6a3037, 
     0x70bc43d9, 0x3b4bb3fa, 0x4a585d0f, 0x58cab8e0, 
     0x2a1f2ff4, 0x59ceade5, 0x228bcdf4, 0x2d0238ee, 
     0x4b30b571, 0x34b8865c, 0x391b17e8, 0x5ff367b5, 
     0x70dbfabc, 0x08d481a1, 0x5462873b, 0x7d4dd4bf,
     0x6a96ceb6, 0x31e29ea8, 0x19d29e1f, 0x7a7d7082, 
     0x7dc1fa60, 0x0eb9819a, 0x11dc28fd, 0x31ba8685, 
     0x5155eb6d, 0x0163fd71, 0x1b4abccf, 0x59adb5e0, 
     0x5b55e0f6, 0x21ccd896, 0x1817e618, 0x4c1224d0, 
     0x5d188c90, 0x62704327, 0x24dcddb0, 0x0737bc84,
     0x3c3ef10c, 0x4768aba4, 0x3439f572, 0x076fa67e, 
     0x7c213200, 0x6d550d5a, 0x67630e33, 0x6cfd2cbd, 
     0x76298efc, 0x3bc5956e, 0x6a4b017c, 0x60c05db2, 
     0x6da83416, 0x041d9f9b, 0x5b3dce34, 0x6b6a2e76, 
     0x12d72135, 0x6d19f731, 0x1d24b4fb, 0x642d0ca2,
     0x6e7df4a3, 0x386f71cb, 0x3ddac282, 0x49d3d599, 
     0x5a3c4a61, 0x55f2a89a, 0x15e5fa69, 0x3754d6f1, 
     0x3862ebc1, 0x3ac2d81a, 0x3e8c9375, 0x74a1dcce, 
     0x022b83be, 0x72c688e8, 0x7c11834c, 0x7e4cb5bf, 
     0x601b9642, 0x6374917f, 0x6b49e27c, 0x5645253e,
     0x1f3a26ee, 0x5594e3f8, 0x370582f0, 0x0ce25b04, 
     0x59b28393, 0x12435124, 0x784c897b, 0x6c89a4c8, 
     0x7f5d4856, 0x15713e76, 0x50b6b16a, 0x6ddb3cf9, 
     0x4de0b041, 0x0e9173ec, 0x37af1292, 0x281cfaa2, 
     0x64841c87, 0x4d950cfc, 0x5f71d193, 0x1ce70848,
     0x0857e516, 0x1dfe6509, 0x1188e516, 0x0a8368d4, 
     0x10c4edf1, 0x0d9a6862, 0x08d01e93, 0x70e08433, 
     0x710ef9e2, 0x741a010f, 0x4725a972, 0x104920d0, 
     0x49aee507, 0x7e2b2c62, 0x1d2b7bd4, 0x2361689a, 
     0x106e7d87, 0x1578054f, 0x0feb0d62, 0x0fcbc5dd,
     0x2ae943c6, 0x60a1becc, 0x7da702d6, 0x78c9f407, 
     0x6f3332b9, 0x35561568, 0x20e6eeaa, 0x53b74f40, 
     0x02eb2264, 0x0058c03d, 0x709e5788, 0x0b43077a, 
     0x1e572546, 0x02273c9f, 0x15c6704f, 0x2f1c1337, 
     0x0fc1a501, 0x1e968ee2, 0x1ffc976b, 0x00d09ee3,
     0x12b08ff2, 0x672240dd, 0x1119bfb3, 0x5c5f74f9, 
     0x654d6d3f, 0x2e453b88, 0x7fc0dd94, 0x75bbeac6, 
     0x43bd40d7, 0x0fabeaf6 };

static int *k_static = m + P;
static int *j_static = m + P-Q;
static int *end          = m + P-1;

static  int        randtbl[ DEG_3 ]              =
  { 0x9a319039, 0x32d9c024, 0x9b663182, 0x5da1f342,
    0xde3b81e0, 0xdf0a6fb5, 0xf103bc02, 0x48f340fb,
    0x7449e56b, 0xbeb1dbb0, 0xab5c5918, 0x946554fd,
    0x8c2e680f, 0xeb3d799f, 0xb11ee0b7, 0x2d436b86,
    0xda672e2a, 0x1588ca88, 0xe369735d, 0x904f35f7,
    0xd7158fd6, 0x6fa6f051, 0x616e6b96, 0xac94efdc,
    0x36413f93, 0xc622c298, 0xf5a42ab8, 0x8a88d77b,
    0xf5ad9d0e, 0x8999220b, 0x27fb47b9 };

static  int                      *fptr    = randtbl + SEP_3;
static  int                      *rptr    = randtbl;
static  int                      *end_ptr = randtbl + DEG_3;
static  int                      *state   = randtbl;
 
static  int                      rand_deg = DEG_3;
static  int                      rand_sep = SEP_3;
 
static int
random_reger( void )

{
            int                      i;
 
            *fptr += *rptr;
            i = (*fptr >> 1) & INT_MAX;      
            /* chucking least random bit */

            if(  ++fptr  >=  end_ptr  )  {
                     fptr = state;   
                     ++rptr;
            }
            else  {
                  if(  ++rptr  >=  end_ptr  )  rptr = state;
            }
            return( i );
}

static void
srandom_reger( unsigned x )

{
            register  int                    i;
                      int random_reger( void );
 
            state[ 0 ] = x;
            for( i = 1; i < rand_deg; i++ )  {
                    state[i] = (1103515245*state[i - 1] + 
                                12345) & INT_MAX;
            }
            fptr = &state[ rand_sep ];
            rptr = &state[ 0 ];
            for( i = 0; i < 10*rand_deg; i++ )  random_reger();
}

void
r250_srandom ( unsigned seed )

{
              register int  i;
 
              srandom_reger( seed );

              for ( i=0; i<P; ++i) m[i] = random_reger();

 k_static = m + P;
 j_static = m + P-Q;
}

int
r250_random( void )
 
{
 
      if (--k_static < m) k_static = end;
      if (--j_static < m) j_static = end;
 
      return (*k_static ^= *j_static);
}
 
float
float_r250( void )
 
{
      if (--k_static < m) k_static = end;
      if (--j_static < m) j_static = end;
 
      return ( (*k_static ^= *j_static) * INVERSE_INT_MAX );
}
 
double
double_r250( void )
 
{
int rnd_int;

      if (--k_static < m) k_static = end;
      if (--j_static < m) j_static = end;
      rnd_int = (*k_static ^= *j_static);
      while (rnd_int==0) {
            if (--k_static < m) k_static = end;
            if (--j_static < m) j_static = end;
            rnd_int = (*k_static ^= *j_static);
      }
      return ( rnd_int * INVERSE_INT_MAX );
}

void
r250_save ( int table[] )

{
      int  i;

      for ( i = 0; i < P; i++ )
 table[i] = m[i];
      
      table[P  ] = k_static - m;
      table[P+1] = j_static - m;
}
 
void
r250_restart ( int table[] )

{
      int  i;

      for ( i = 0; i < P; i++ )
 m[i] = table[i];
      
      k_static = m + table[P  ];
      j_static = m + table[P+1];
}
 
double
gauss_r250 ( double mean, double variance )
 
{
      static int           empty = 1;
      static double        store;
      double               v1, v2, r, root;

      if ( empty ) {
 
 do { /* Inline double_random() */

      if (--k_static < m) k_static = end;
      if (--j_static < m) j_static = end;
      v1 = (*k_static ^= *j_static) * INVERSE_INT_MAX;

      v1 = 2.0 * v1 - 1.0;

      if (--k_static < m) k_static = end;
      if (--j_static < m) j_static = end;
      v2 = (*k_static ^= *j_static) * INVERSE_INT_MAX;

      v2 = 2.0 * v2 - 1.0;

      r  = v1 * v1 + v2 * v2;

 } while ( r > 1.0 );

 root  = variance * sqrt( -2.0 * log(r) / r );
 store = v2 * root + mean;
 empty = 0;

 return ( v1 * root + mean );

      } else {

 empty = 1;
 return ( store );

      }
 
}



# ifndef PROBLEM_R250
      # define MAX(a,b) ((a)>(b)?(a):(b))
# else
      # define MAX(a,b) ((((a)>(b))&&((a)<(int *)\
0x80000000  ))?(a):(b))
# endif

void
r250_vector ( register int *x, int n)
 
{
      register int         *j, *k, *l;
 
      j = j_static;
      k = k_static;

      while ( n > 0 ) {

 if ( k == m ) k = m + P;
 l = MAX( k-n, m+Q );
 if ( k > l ) {
      n -= k-l;
#pragma ivdep
      while ( k > l )
         *x++ = (*--k ^= *--j);
 }

 if ( j == m ) j = m + P;
 l = MAX( j-n, m+P-Q );
 if ( j > l ) {
      n -= j-l;
#pragma ivdep
      while ( j > l )
         *x++ = (*--k ^= *--j);
 }

      }

      j_static = j;
      k_static = k;

}



void
float_r250_vector ( register float *x, int n)
 
{
      register int         *j, *k, *l;
 
      j = j_static;
      k = k_static;

      while ( n > 0 ) {

 if ( k == m ) k = m + P;
 l = MAX( k-n, m+Q );
 if ( k > l ) {
      n -= k-l;
#pragma ivdep
      while ( k > l )
         *x++ = (*--k ^= *--j) * INVERSE_INT_MAX;
 }

 if ( j == m ) j = m + P;
 l = MAX( j-n, m+P-Q );
 if ( j > l ) {
      n -= j-l;
#pragma ivdep
      while ( j > l )
         *x++ = (*--k ^= *--j) * INVERSE_INT_MAX;
 }

      }

      j_static = j;
      k_static = k;

}

#ifndef SERIOUS_R250
void
double_r250_vector ( register double *x, int n)
 
{
      register int         *j, *k, *l;
      int                                rnd; 
 
      j = j_static;
      k = k_static;

      while ( n > 0 ) {

 if ( k == m ) k = m + P;
 l = MAX( k-n, m+Q );
 if ( k > l ) {
      n -= k-l;
#pragma ivdep
      while ( k > l ) {
         rnd = (*--k ^= *--j);
         if (rnd==0) rnd=1;
         *x++ = rnd * INVERSE_INT_MAX;
      }
 }

 if ( j == m ) j = m + P;
 l = MAX( j-n, m+P-Q );
 if ( j > l ) {
      n -= j-l;
#pragma ivdep
      while ( j > l ) {
         rnd = (*--k ^= *--j);
         if (rnd==0) rnd=1;
         *x++ = rnd * INVERSE_INT_MAX;
      }
 }

      }

      j_static = j;
      k_static = k;

}
#else
void
double_r250_vector ( register double *x, int n)
{
double *ptr,rnd;
int      i,rnd_int;

    for (i=0,ptr=x;i<n;i++) {
      if (--k_static < m) k_static = end;
      if (--j_static < m) j_static = end;
      rnd_int = (*k_static ^= *j_static);
      while (rnd_int==0) {
            if (--k_static < m) k_static = end;
            if (--j_static < m) j_static = end;
            rnd_int = (*k_static ^= *j_static);
      }
      rnd  = rnd_int * INVERSE_INT_MAX;
      *x++ = rnd;
    }  
}
#endif



/***************************************************/
/*                                                 */
/* File: r250.c                                    */
/*                                                 */
/* A simple multiplicative generator from          */
/* Kalos & Whitlock.                               */
/*                                                 */
/***************************************************/

static int       rand_seed = 1;

int
int_rand( void )

{
rand_seed = ( 1812433253 * rand_seed + 314159265 ) & INT_MAX;
return  ( rand_seed );
}

float
float_rand( void )

{
rand_seed = ( 1812433253 * rand_seed + 314159265 ) & INT_MAX;
return  ( rand_seed * INVERSE_INT_MAX );
}


void
save_rand( int seed )

{
      rand_seed = seed;
}

int
rand_save( void )

{
      return ( rand_seed );
}
