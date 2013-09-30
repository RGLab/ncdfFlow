#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include "hdfFlow.h"
#define ERRCODE 2

//#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define ERR(e) {REprintf("hdf Error: %s\n", nc_strerror(e));SEXP k = allocVector(LGLSXP,1); LOGICAL(k)[0]=FALSE;return(k);}


#define CACHE_SIZE 64000000
#define CACHE_NELEMS 1009
#define CACHE_PREEMPTION .75
#define DEFLATE_LEVEL 2

SEXP createFile(SEXP _fileName, SEXP _X, SEXP _Y, SEXP _Z, SEXP _metaSize,SEXP _compress) {
    //Rprintf("createFile\n");
    int NDIMS = 3 ; 
    int dimids[NDIMS];
    int X = INTEGER(_X)[0];
    int Y = INTEGER(_Y)[0];
    int Z = INTEGER(_Z)[0];
    //dimension size for metadata is as double as the initial input size
    //for the future possible increase of metaData
    int metaSize=INTEGER(_metaSize)[0]*2;
    int M=1;//meta
    int ncid, retval, varid, x_dimid, y_dimid, z_dimid;
    int m_dimid,m_varid;//metadata
    int compress = LOGICAL(_compress)[0];
    SEXP k = allocVector(LGLSXP,1);
    size_t chunksize[] = {1, 1, X};
    
    if ((retval = nc_create( translateChar(STRING_ELT(_fileName, 0)), NC_NETCDF4, &ncid)))
      ERR(retval);


    if ((retval = nc_def_dim(ncid, "meta",metaSize, &m_dimid)))
        ERR(retval);

    if ((retval = nc_def_dim(ncid, "event", X, &x_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "channel", Y, &y_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "sample", Z , &z_dimid))) // was NC_UNLIMITED Z 
        ERR(retval);

    dimids[0] = z_dimid;
    dimids[1] = y_dimid;
    dimids[2] = x_dimid;

    if ((retval = nc_def_var(ncid, "metaData", NC_BYTE, M, &m_dimid, &m_varid)))//metadata
        ERR(retval);

    if ((retval = nc_def_var(ncid, "exprsMat", NC_DOUBLE, NDIMS, dimids, &varid)))
        ERR(retval);
    
    if (( retval = nc_def_var_chunking(ncid, varid, NC_CHUNKED, chunksize)))
        ERR(retval);
    	/*use default cache settings since we are not using it at the moment*/
//    if (( retval = nc_set_var_chunk_cache(ncid, varid, CACHE_SIZE, CACHE_NELEMS, CACHE_PREEMPTION)))
//        ERR(retval);
    
    if(compress) {
        if (( retval = nc_def_var_deflate(ncid, varid, 0, 1, DEFLATE_LEVEL)))
            ERR(retval);
    }

    if ((retval = nc_enddef(ncid)))
        ERR(retval);

    //Save variable attributes 
    //No of channels stored
    if((retval = nc_put_att_int(ncid, varid, "channelCount", NC_INT, 1, &Y)))
        ERR(retval);
    // No of samples 
    if((retval = nc_put_att_int(ncid, varid, "sampleCount", NC_INT, 1, &Z)))
        ERR(retval);
 
    int *eCount = (int *) R_alloc( sizeof(int), Z);
    for(int i =0 ; i < Z; i++) {
        eCount[i] = 0;
    }
    //No of events stored
    if((retval = nc_put_att_int(ncid, varid, "eventCount", NC_INT, Z, eCount)))
        ERR(retval);
    //attribute to record the actual size of metaData variable
    metaSize=0;
    if((retval = nc_put_att_int(ncid, m_varid, "metaSize", NC_INT, 1, &metaSize)))
            ERR(retval);

    //attribute to record the maximum size of metaData variable
//     if((retval = nc_put_att_int(ncid, m_varid, "maxMetaSize", NC_INT, 1, &metaSize)))
//                ERR(retval);

    if ((retval = nc_close(ncid)))
        ERR(retval);
    LOGICAL(k)[0]=TRUE;
    return(k);
}


/*each channel  is stored as a data chunk(size= events), which is
indexed in file and fast to access as a whole unit(slice),so the C API for accessing
events  are not provided due to the inefficient IO.   */
SEXP writeSlice(SEXP _fileName, SEXP _mat, SEXP _chIndx, SEXP _sample) {
    //Rprintf("writeSlice\n");
    int retval, ncid, varid, nRow, nCol, sample;
    SEXP Rdim, k = allocVector(LGLSXP,1);
    Rdim = getAttrib(_mat, R_DimSymbol);
    nRow = INTEGER(Rdim)[0];
    nCol = INTEGER(Rdim)[1];
    sample = INTEGER(_sample)[0]-1;  // R to C indexing


    double *mat = REAL(_mat);
    int * chIndx = INTEGER(_chIndx);
	int chCount = length(_chIndx);
    

    if ((retval = nc_open(translateChar(STRING_ELT(_fileName, 0)), NC_WRITE, &ncid)))
        ERR(retval);
    

    if((retval = nc_inq_varid (ncid, "exprsMat", &varid)))
        ERR(retval);

    for(unsigned i=0;i<chCount;i++)
	{
		int colStart = chIndx[i]-1;
		size_t start[] = {sample, colStart, 0};
		size_t count[] = {1, 1, nRow};
		double * vec = mat+i*nRow;
		if((retval = nc_put_vara_double(ncid, varid, start, count, vec)))
				ERR(retval);

      }

    int sampCount;

    if((retval = nc_get_att_int(ncid, varid, "sampleCount", &sampCount)))
        ERR(retval);

    int len = sampCount;
    if(sample >= sampCount) {
        len = sample;
    }

    int *eCount = (int *) R_alloc( sizeof(int), len);
    //Rprintf("Get eventCount attribute\n");
    if((retval = nc_get_att_int(ncid, varid, "eventCount", eCount)))
        ERR(retval);
    eCount[sample] = nRow;


    //Rprintf("Redefine ncid\n");
	if((retval = nc_redef(ncid)))
		ERR(retval);
    if((retval = nc_put_att_int(ncid, varid, "eventCount", NC_INT, len, eCount)))
		ERR(retval);
    //Rprintf("Finished defining ncid\n");
    if((retval = nc_enddef(ncid)))
        ERR(retval);
    //nc_sync(ncid);
   // Rprintf("Closing Ncfile\n");
    if ((retval = nc_close(ncid)))
        ERR(retval);
    LOGICAL(k)[0]=TRUE;
    return(k);

}

SEXP readSlice(SEXP _fileName, SEXP _chIndx, SEXP _sample ) {

    int retval, ncid, varid, nRow;
    SEXP ans, dnms;
    int sample = INTEGER(_sample)[0]-1;  // R to C indexing
    int * chIndx = INTEGER(_chIndx);
    int chCount = length(_chIndx);


    if ((retval = nc_open(translateChar(STRING_ELT(_fileName, 0)), NC_NOWRITE,&ncid))){

        ERR(retval);
    }
        
    

    if((retval = nc_inq_varid (ncid, "exprsMat", &varid)))
        ERR(retval);  

    int sampCount;

    if((retval = nc_get_att_int(ncid, varid, "sampleCount", &sampCount)))
        ERR(retval);
 
    int *eCount = (int *) R_alloc( sizeof(int), sampCount);

    if((retval = nc_get_att_int(ncid, varid, "eventCount", eCount)))
        ERR(retval);
    nRow = eCount[sample] ;
    
    PROTECT(ans = allocVector(REALSXP, nRow*chCount));
    double *mat = REAL(ans);
    for(unsigned i=0;i<chCount;i++){
    	int colStart = chIndx[i]-1;
    	size_t start[] = {sample, colStart, 0};
		size_t count[] = {1, 1, nRow};
		double * vec = mat+i*nRow;
		if((retval = nc_get_vara_double(ncid, varid, start, count, vec)))
			ERR(retval);
    }

    if ((retval = nc_close(ncid)))
        ERR(retval);
    PROTECT(dnms = allocVector(INTSXP, 2));
    INTEGER(dnms)[0] = nRow;
    INTEGER(dnms)[1]=  chCount;
    setAttrib(ans,R_DimSymbol, dnms);
    UNPROTECT(2);
    return(ans);
}



