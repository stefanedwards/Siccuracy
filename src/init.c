#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    subroutine

  Most likely possible values need to be added below.
*/


/* .Fortran calls */
extern void F77_NAME(cbindsnpsrwrapper)(int *, char *, char *, int *, int *, int *, int *, int *, int *, int *, char *, int *);
extern void F77_NAME(convert_phase)(char *, char *, int *, int *, int *, int *, int *, char *);
extern void F77_NAME(convertplinka)(char *, char *, int *, int *, int *, int *, int *, int *);
extern void F77_NAME(convertplinkrwrapper)(char *, int *, int *, int *, char *, char *, int *, int *, int *, int *, int *, double *, int *, int *, int *);
extern void F77_NAME(get_nlines)(char *, int *, int *);
extern void F77_NAME(heterozygosity)(char *, int *, int *, int *, int *, double *, double *, double *, int *);
extern void F77_NAME(imp_acc)(char *, char *, int *, int *, int *, int *, double *, double *, int *, double *, double *, double *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void F77_NAME(imp_acc_fast)(char *, char *, int *, int *, int *, int *, double *, double *, int *, double *, double *, double *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void F77_NAME(masksnps2)(char *, char *, int *, int *, int *, char *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void F77_NAME(rbindsnps)(char *, char *, char *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, char *, int *, int *, int *);
extern void F77_NAME(readplinksimple)(char *, char *, int *, int *, int *, int *, int *, double *, int *, int *, int *);

static const R_FortranMethodDef FortranEntries[] = {
    {"cbindsnpsrwrapper",    (DL_FUNC) &F77_NAME(cbindsnpsrwrapper),    12},
    {"convert_phase",        (DL_FUNC) &F77_NAME(convert_phase),         8},
    {"convertplinka",        (DL_FUNC) &F77_NAME(convertplinka),         8},
    {"convertplinkrwrapper", (DL_FUNC) &F77_NAME(convertplinkrwrapper), 15},
    {"get_nlines",           (DL_FUNC) &F77_NAME(get_nlines),            3},
    {"heterozygosity",       (DL_FUNC) &F77_NAME(heterozygosity),        9},
    {"imp_acc",              (DL_FUNC) &F77_NAME(imp_acc),              24},
    {"imp_acc_fast",         (DL_FUNC) &F77_NAME(imp_acc_fast),         24},    
    {"masksnps2",            (DL_FUNC) &F77_NAME(masksnps2),            18},
    {"rbindsnps",            (DL_FUNC) &F77_NAME(rbindsnps),            18},
    {"readplinksimple",      (DL_FUNC) &F77_NAME(readplinksimple),      11},
    {NULL, NULL, 0}
};

void R_init_Siccuracy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
