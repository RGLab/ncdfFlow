#include <Rconfig.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "../inst/include/ncdfFlowAPI.h"


void R_init_ncdfFlow(DllInfo *dll)
{

	R_RegisterCCallable("ncdfFlow","ncdfFlow_open_hdf",(DL_FUNC) &_ncdfFlow_open_hdf);
	R_RegisterCCallable("ncdfFlow","ncdfFlow_get_event_number",(DL_FUNC) &_ncdfFlow_get_event_number);
	R_RegisterCCallable("ncdfFlow","ncdfFlow_readSlice",(DL_FUNC) &_ncdfFlow_readSlice);
	R_RegisterCCallable("ncdfFlow","ncdfFlow_close_hdf",(DL_FUNC) &_ncdfFlow_close_hdf);
}
