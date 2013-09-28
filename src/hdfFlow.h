#include "hdf5.h"

herr_t _createFile(const char * fName, unsigned nSample, unsigned nChnl, unsigned nEvt);
void _writeSlice(const char * fName, double * mat, unsigned chnlIndx,  unsigned sampleIndx);
double * _readSlice(const char * fName, unsigned chnlIndx, unsigned sampleIndx);
