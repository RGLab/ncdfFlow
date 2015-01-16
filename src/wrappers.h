#ifndef WRAPPERS_H_
#define WRAPPERS_H_

#include "hdf5.h"
#include <stdlib.h>
#include <stdio.h>
#define DATASETNAME3d "/exprsMat"
#define MAXLEN 100//max character length of sample index

#define TRUE            1
#define FALSE           0

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>





herr_t custom_print_cb(hid_t estack, void *client_data);
herr_t my_hdf5_error_handler(unsigned n, const H5E_error2_t *err_desc, void *client_data);
SEXP createFile(SEXP _fileName, SEXP _nEvent, SEXP _nChannel, SEXP _nSample, SEXP _dim, SEXP _ratio);

SEXP writeSlice(SEXP _fileName, SEXP _mat, SEXP _chIndx, SEXP _sampleIndx, SEXP _ratio);

SEXP readSlice(SEXP _fileName, SEXP _chIndx, SEXP _sampleIndx, SEXP _colnames);


#endif /* WRAPPERS_H_ */
