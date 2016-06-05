#ifndef __PYX_HAVE__det
#define __PYX_HAVE__det


#ifndef __PYX_HAVE_API__det

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

__PYX_EXTERN_C DL_IMPORT(double) det(double **, int);

#endif /* !__PYX_HAVE_API__det */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initdet(void);
#else
PyMODINIT_FUNC PyInit_det(void);
#endif

#endif /* !__PYX_HAVE__det */
