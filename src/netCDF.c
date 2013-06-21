#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <hdf5.h>
#define ERRCODE 2

//#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define ERR(e) {REprintf("netCDF Error: %s\n", nc_strerror(e));SEXP k = allocVector(LGLSXP,1); LOGICAL(k)[0]=FALSE;return(k);}


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

    if ((retval = nc_def_var(ncid, "exprsMat", NC_FLOAT, NDIMS, dimids, &varid)))
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


/*each sample is stored as a data chunk(size=channel x events), which is
indexed in file and fast to access as a whole unit(slice),so the C API for accessing
events or channels are not provided due to the inefficient IO.   */
SEXP writeSlice(SEXP _fileName, SEXP _mat, SEXP _sample) {
    //Rprintf("writeSlice\n");
    int retval, ncid, varid, nRow, nCol, sample;
    SEXP Rdim, k = allocVector(LGLSXP,1);
    Rdim = getAttrib(_mat, R_DimSymbol);
    nRow = INTEGER(Rdim)[0];
    nCol = INTEGER(Rdim)[1];
    sample = INTEGER(_sample)[0]-1;  // R to C indexing
    size_t start[] = {sample, 0, 0};
    size_t count[] = {1, nCol, nRow};
    double *mat = REAL(_mat);
    
    //Rprintf("Opening file\n");
    if ((retval = nc_open(translateChar(STRING_ELT(_fileName, 0)), NC_WRITE, &ncid)))
        ERR(retval);
    
    //Rprintf("Retrieving variable\n");
    if((retval = nc_inq_varid (ncid, "exprsMat", &varid)))
        ERR(retval);
    // check if max rows of ncdf file is exceeded at R level
    //Rprintf("Write variable\n");
    if((retval = nc_put_vara_double(ncid, varid, start, count, mat)))
        ERR(retval);

    int sampCount;
    //Rprintf("Get sampleCount attribute\n");
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
    //Rprintf("Writing eventCount\n");
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

SEXP readEventCounts(SEXP _fileName) {
    //Rprintf("readEventCounts\n");
    int retval, ncid, varid;
    SEXP ans;
    if ((retval = nc_open(translateChar(STRING_ELT(_fileName, 0)), NC_NOWRITE, &ncid)))
        ERR(retval);

    if((retval = nc_inq_varid (ncid, "exprsMat", &varid)))
        ERR(retval);  
     
    int sampCount;
    if((retval = nc_get_att_int(ncid, varid, "sampleCount", &sampCount)))
        ERR(retval);
    
    PROTECT(ans = allocVector(INTSXP, sampCount));
    int *eCount = INTEGER(ans);
    if((retval = nc_get_att_int(ncid, varid, "eventCount", eCount)))
        ERR(retval);

    if ((retval = nc_close(ncid)))
        ERR(retval);
    UNPROTECT(1);
    return(ans);
}


SEXP writeMeta(SEXP _fileName, SEXP _meta,SEXP _start,SEXP _count) {
    //Rprintf("writeMeta\n");
	SEXP k = allocVector(LGLSXP,1);
    int retval, ncid, varid;
//    int dimid;
    size_t start =INTEGER(_start)[0]-1 ;
    size_t count = INTEGER(_count)[0];
//    int maxMetaSize;
    unsigned char * meta =RAW(_meta);

    if ((retval = nc_open(translateChar(STRING_ELT(_fileName, 0)), NC_WRITE, &ncid)))
        ERR(retval);

    if((retval = nc_inq_varid (ncid, "metaData", &varid)))
        ERR(retval);
//     check if max rows of ncdf file is exceeded at R level
    if((retval = nc_put_vara(ncid, varid, &start,&count,meta)))
        ERR(retval);

    int metaSize=count;

//    if((retval = nc_inq_dimid (ncid, "meta", &dimid)))
//        ERR(retval);
//    if((retval = nc_inq_dimlen (ncid,dimid, &maxMetaSize)))
//          ERR(retval);
//    if((retval = nc_get_att_int(ncid, varid, "maxMetaSize", &maxMetaSize)))
//            ERR(retval);
//    if(metaSize>maxMetaSize)
//    {
//    	printf("Error: %s\n", "Meta data to write is too big!");
//    	k=-1;
//    	return(k);
//    }

    if((retval = nc_redef(ncid)))
		ERR(retval);
    //update the actual size of metaData variable
	if((retval = nc_put_att_int(ncid, varid, "metaSize", NC_INT, 1, &metaSize)))
		ERR(retval);
	if((retval = nc_enddef(ncid)))
		ERR(retval);

    if ((retval = nc_close(ncid)))
        ERR(retval);
    LOGICAL(k)[0]=TRUE;
    return(k);

}

SEXP readMeta(SEXP _fileName) {

	SEXP ans;
    int retval, ncid, varid,metaSize;
    size_t count,start =0;


    if ((retval = nc_open(translateChar(STRING_ELT(_fileName, 0)), NC_NOWRITE, &ncid)))
        ERR(retval);

    if((retval = nc_inq_varid (ncid, "metaData", &varid)))
        ERR(retval);

    if((retval = nc_get_att_int(ncid, varid, "metaSize",&metaSize)))
    		ERR(retval);

    PROTECT(ans = allocVector(RAWSXP,metaSize));
    unsigned char * meta =RAW(ans);

    count=metaSize;
    //     check if max rows of ncdf file is exceeded at R level
    if((retval = nc_get_vara(ncid, varid, &start,&count,meta)))
        ERR(retval);

    if ((retval = nc_close(ncid)))
        ERR(retval);
    UNPROTECT(1);

    return(ans);

}


SEXP createIndiceFile(SEXP _fileName, SEXP _eventCount, SEXP _totalNodeCount) {
    SEXP k = allocVector(LGLSXP,1);
	int NDIMS = 2 ;
	int dimids[NDIMS];
	/*
	 * nBytes count is the actual storage length,which can have extra bits that should
	 * not be converted to logical values,so eventCount is necessary information for bit vector routine
	 *  to retrieve exact length of logical vector
	 */
	int eventCount = INTEGER(_eventCount)[0];
    int sizeInBytes=ceil((float)eventCount/8);
    int totalNodeCount = INTEGER(_totalNodeCount)[0];
    int ncid, retval, varid,x_dimid, y_dimid;
    size_t chunksize[] = {1,sizeInBytes};

    if ((retval = nc_create( translateChar(STRING_ELT(_fileName, 0)), NC_NETCDF4, &ncid)))
      ERR(retval);

    if ((retval = nc_def_dim(ncid, "indice", sizeInBytes, &x_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "Node", totalNodeCount, &y_dimid)))
        ERR(retval);


    dimids[0] = y_dimid;
    dimids[1] = x_dimid;

    if ((retval = nc_def_var(ncid, "IndiceMat", NC_BYTE, NDIMS, dimids, &varid)))
        ERR(retval);

    if (( retval = nc_def_var_chunking(ncid, varid, NC_CHUNKED, chunksize)))
        ERR(retval);

    if (( retval = nc_set_var_chunk_cache(ncid, varid, CACHE_SIZE, CACHE_NELEMS, CACHE_PREEMPTION)))
        ERR(retval);

    //compress
    if (( retval = nc_def_var_deflate(ncid, varid, 0, 1, DEFLATE_LEVEL)))
            ERR(retval);

    //length of actual bits that are used ,equal to the number of events
    if((retval = nc_put_att_int(ncid, varid, "bitlen", NC_INT, 1, &eventCount)))
            ERR(retval);
    // No of gates/nodes
    if((retval = nc_put_att_int(ncid, varid, "totalNodeCount", NC_INT, 1, &totalNodeCount)))
            ERR(retval);
    // No of gates/nodes
    if((retval = nc_put_att_int(ncid, varid, "sizeInBytes", NC_INT, 1, &sizeInBytes)))
               ERR(retval);


    //nbitset store the bits that are 1,which would be used in bitvector c routine
    int *nbitset = (int *) R_alloc( sizeof(int), totalNodeCount);
    for(int i =0 ; i < totalNodeCount; i++) {
    	nbitset[i] = 0;
    }
    //No of bits stored
    if((retval = nc_put_att_int(ncid, varid, "nbitset", NC_INT,totalNodeCount, nbitset)))
	  ERR(retval);

    if ((retval = nc_close(ncid)))
        ERR(retval);

    LOGICAL(k)[0]=TRUE;

    return(k);


}

SEXP writeIndice(SEXP _fileName, SEXP _IndiceMat, SEXP _NodeIndStart) {

    int retval;
    //return logical value
    SEXP k = allocVector(LGLSXP,1);
    //get dimention info from input matrix
    SEXP Rdim = getAttrib(_IndiceMat, R_DimSymbol);
    int nLength = INTEGER(Rdim)[0];
    int nNodes = INTEGER(Rdim)[1];//number of nodes to write

    //get start position from input
    int nodeIndStart = INTEGER(_NodeIndStart)[0]-1;  // R to C indexing
    size_t start[] = {nodeIndStart, 0};
    size_t count[] = {nNodes,nLength};

    //open file
	int ncid;
	if ((retval = nc_open(translateChar(STRING_ELT(_fileName, 0)), NC_WRITE, &ncid)))
		ERR(retval);
	//get variable id
	int varid;
	if((retval = nc_inq_varid (ncid, "IndiceMat", &varid)))
		ERR(retval);

    //retrieve the length of indice vector
    int sizeInBytes;
    if((retval = nc_get_att_int(ncid, varid, "sizeInBytes", &sizeInBytes)))
              ERR(retval);
    //check the validity of input matrix size
    if(sizeInBytes!=nLength)
    {
    	REprintf("Error: %s\n", "logical vector size is not the same as the one in cdf !");
        LOGICAL(k)[0]=FALSE;//set return as true
        return(k);
    }



    //write matrix
    unsigned char * Indice =RAW(_IndiceMat);
    if((retval = nc_put_vara(ncid, varid, start, count, Indice)))
        ERR(retval);

    //retrieve total number of nodes
    int totalNodeCount;
    if((retval = nc_get_att_int(ncid, varid, "totalNodeCount", &totalNodeCount)))
         ERR(retval);


    //get attribute info from input matrix
    SEXP Rattr= getAttrib(_IndiceMat, install("nbitset"));
    int *nbitset_input=INTEGER(Rattr);
    //retrieve the nbitset array
    int *nbitset = (int *) R_alloc( sizeof(int), totalNodeCount);
    if((retval = nc_get_att_int(ncid, varid, "nbitset", nbitset)))
         ERR(retval);
    //update part of it
    for(int i =0 ; i < nNodes; i++) {
       	nbitset[i+nodeIndStart] = nbitset_input[i];
       }
    //write the whole array back
    if((retval = nc_redef(ncid)))
         ERR(retval);
    if((retval = nc_put_att_int(ncid, varid, "nbitset", NC_INT, totalNodeCount, nbitset)))
    		ERR(retval);
    if((retval = nc_enddef(ncid)))
         ERR(retval);

    //close file
    if ((retval = nc_close(ncid)))
        ERR(retval);

    LOGICAL(k)[0]=TRUE;//set return as true
    return(k);

}

SEXP readIndice(SEXP _fileName,SEXP _NodeIndStart,SEXP _nNode) {

    int retval;
//    SEXP k = allocVector(LGLSXP,1);

    int nodeIndStart = INTEGER(_NodeIndStart)[0]-1;  // R to C indexing
    int nNode = INTEGER(_nNode)[0];  // R to C indexing

    //open file
    int ncid;
    if ((retval = nc_open(translateChar(STRING_ELT(_fileName, 0)), NC_NOWRITE, &ncid)))
        ERR(retval);
    //get variable ID
    int varid;
    if((retval = nc_inq_varid (ncid, "IndiceMat", &varid)))
        ERR(retval);
    //get variable length
    int sizeInBytes;
    if((retval = nc_get_att_int(ncid, varid, "sizeInBytes", &sizeInBytes)))
                  ERR(retval);

    //retrieve indices in raw bytes and store them in matrix
    const size_t start[] = {nodeIndStart, 0};
    const size_t count[] = {nNode,sizeInBytes};
    SEXP ans;
    PROTECT(ans = allocVector(RAWSXP,sizeInBytes*nNode));
    unsigned char * mat =RAW(ans);
    if((retval = nc_get_vara(ncid, varid, start,count,mat)))
            ERR(retval);

    //construct dimentions and add it to matrix
    SEXP dnms;
    PROTECT(dnms = allocVector(INTSXP, 2));
    INTEGER(dnms)[0] = sizeInBytes;
    INTEGER(dnms)[1]=nNode ;
    setAttrib(ans,R_DimSymbol, dnms);

    //construct attribute bitlen and add it to matrix
    int bitlen;
    if((retval = nc_get_att_int(ncid, varid, "bitlen",&bitlen)))
         ERR(retval);
    SEXP btln;
    PROTECT(btln = allocVector(INTSXP, 1));
    INTEGER(btln)[0]=bitlen;
    setAttrib(ans, install("bitlen"), btln);

    //construct attribute nbitset and add it to matrix
    int *nbitset = (int *) R_alloc( sizeof(int), sizeInBytes);
    if((retval = nc_get_att_int(ncid, varid, "nbitset",nbitset)))
     	ERR(retval);
    SEXP btcnt;
    PROTECT(btcnt = allocVector(INTSXP, nNode));
    for(int i=0;i<nNode;i++)
    	INTEGER(btcnt)[i]=nbitset[i+nodeIndStart];
    setAttrib(ans, install("nbitset"), btcnt);

    //close file
    if ((retval = nc_close(ncid)))
            ERR(retval);

    UNPROTECT(4);

    return(ans);


}
