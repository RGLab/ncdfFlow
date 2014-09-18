#include "wrappers.h"

static const R_CallMethodDef CallEntries[] = {
    {"createFile", (DL_FUNC)&createFile, 6},
    {"writeSlice", (DL_FUNC)&writeSlice, 5},
    {"readSlice", (DL_FUNC)&readSlice, 4},
    {NULL, NULL, 0}
};

void R_init_ncdfFlow(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);

}
