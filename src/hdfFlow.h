#ifndef HDFFLOW_H_
#define HDFFLOW_H_

#include "hdf5.h"
#include <RcppArmadillo.h>

#define DATASETNAME3d "/exprsMat"

#define TRUE            1
#define FALSE           0

#define MSG_SIZE       1024

#include <vector>
#include <string>

herr_t custom_print_cb(hid_t estack, void *client_data);

herr_t my_hdf5_error_handler(unsigned n, const H5E_error2_t *err_desc, void *client_data);

bool createFile(std::string fileName, int nEvent, int nChannel, int nSample, int nDim, int nCompressionRatio);

bool writeSlice(std::string filename, Rcpp::NumericMatrix _mat, std::vector<int> chIndx, int sampleIndx, int nRatio);

void open_hdf(std::string filename, unsigned flags, hid_t & fileid, hid_t & dataset, hid_t & dataspace, bool & is3d);

unsigned get_event_number(hid_t fileid, hid_t & dataset, hid_t & dataspace, unsigned sampleIndx, bool is3d);

void readSlice_cpp(hid_t fileid, hid_t dataset, hid_t dataspace, std::vector<unsigned> chIndx, unsigned sampleIndx, unsigned nEvents, double * data_out, bool is3d);

void readSlice(hid_t fileid, hid_t dataset, hid_t dataspace, std::vector<unsigned> chIndx, unsigned sampleIndx, unsigned nEvents, Rcpp::NumericVector data_out, bool is3d);

void close_hdf(hid_t fileid);

#endif /* HDFFLOW_H_ */
