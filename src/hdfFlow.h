#include "hdf5.h"
#include <R.h>
#include <stdlib.h>
#define DATASETNAME3d "/exprsMat"


#define TRUE            1
#define FALSE           0

herr_t _createFile3d(const char * fName, unsigned nSample, unsigned nChnl, unsigned nEvt);
herr_t _writeSlice(const char * fName, double * mat, unsigned nEvents, unsigned * chnlIndx, unsigned chCount,  unsigned sampleIndx);
herr_t _readSlice(const char * fName, unsigned * chnlIndx, unsigned chCount, unsigned sampleIndx, double * data_out);


herr_t _createFile2d(const char * fName);
