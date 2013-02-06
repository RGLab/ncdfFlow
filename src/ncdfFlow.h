
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
SEXP ncdf_bitarray_set(SEXP bits, SEXP _indx);
SEXP ncdf_bitarray_getSetBitPos(SEXP bits);
SEXP ncdf_bitarray_getBitStatus(SEXP bits);
SEXP ncdf_bitarray_Flip(SEXP bits);

SEXP createFile(SEXP _fileName, SEXP _X, SEXP _Y, SEXP _Z, SEXP _metaSize,SEXP _compress);
SEXP writeSlice(SEXP _fileName, SEXP _mat, SEXP _sample);
SEXP readSlice(SEXP _fileName, SEXP _y, SEXP _sample );
SEXP writeMeta(SEXP _fileName, SEXP _meta,SEXP _start,SEXP _count);
SEXP readMeta(SEXP _fileName);
SEXP createIndiceFile(SEXP _fileName, SEXP _eventCount, SEXP _totalNodeCount);
SEXP writeIndice(SEXP _fileName, SEXP _IndiceMat, SEXP _NodeIndStart);
SEXP readIndice(SEXP _fileName,SEXP _NodeIndStart,SEXP _nNode);
