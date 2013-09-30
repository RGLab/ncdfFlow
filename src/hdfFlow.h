#include "hdf5.h"
#include <stdlib.h>
#define DATASETNAME "/exprsMat"

herr_t _createFile(const char * fName, unsigned nSample, unsigned nChnl, unsigned nEvt);
herr_t _writeSlice(const char * fName, double * mat, unsigned nEvents, unsigned * chnlIndx, unsigned chCount,  unsigned sampleIndx);
herr_t _readSlice(const char * fName, unsigned * chnlIndx, unsigned chCount, unsigned sampleIndx, double * data_out);
