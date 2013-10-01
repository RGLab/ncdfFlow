#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include "hdfFlow.h"

#define ERR(e) {REprintf("hdf Error: %s\n", HEstring(e));SEXP k = allocVector(LGLSXP,1); LOGICAL(k)[0]=FALSE;return(k);}


SEXP createFile(SEXP _fileName, SEXP _X, SEXP _Y, SEXP _Z, SEXP _metaSize,SEXP _compress) {

    int X = INTEGER(_X)[0];
    int Y = INTEGER(_Y)[0];
    int Z = INTEGER(_Z)[0];
    int compress = LOGICAL(_compress)[0];

    int retval = _createFile(translateChar(STRING_ELT(_fileName, 0)), X, Y, Z);
    if(retval < 0)
      ERR(retval);

    return(R_NilValue);
}


SEXP writeSlice(SEXP _fileName, SEXP _mat, SEXP _chIndx, SEXP _sample) {

    int retval, nRow, nCol, sample;
    SEXP Rdim = getAttrib(_mat, R_DimSymbol);
    nRow = INTEGER(Rdim)[0];
    nCol = INTEGER(Rdim)[1];
    sample = INTEGER(_sample)[0];


    double *mat = REAL(_mat);
    int * chIndx = INTEGER(_chIndx);
	int chCount = length(_chIndx);
    
	retval = _writeSlice(translateChar(STRING_ELT(_fileName, 0)), mat, nRow, chIndx, chCount, sample);
	if(retval < 0)
        ERR(retval);

    return(R_NilValue);

}

SEXP readSlice(SEXP _fileName, SEXP _chIndx, SEXP _sample ) {

    int retval, ncid, varid, nRow;
    SEXP ans, dnms;
    int sample = INTEGER(_sample)[0];
    int * chIndx = INTEGER(_chIndx);
    int chCount = length(_chIndx);

    PROTECT(ans = allocVector(REALSXP, nRow*chCount));
    double *mat = REAL(ans);

    retval = _readSlice(translateChar(STRING_ELT(_fileName, 0)), chIndx, chCount, sample, mat);
	if(retval < 0)
    	ERR(retval);
    

    PROTECT(dnms = allocVector(INTSXP, 2));
    INTEGER(dnms)[0] = nRow;
    INTEGER(dnms)[1]=  chCount;
    setAttrib(ans,R_DimSymbol, dnms);
    UNPROTECT(2);
    return(ans);
}



