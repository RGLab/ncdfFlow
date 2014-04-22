#include "wrappers.h"

static const R_CallMethodDef CallEntries[] = {
    {"bitarray_set", (DL_FUNC)&bitarray_set, 2},
    {"bitarray_getSetBitPos", (DL_FUNC)&bitarray_getSetBitPos, 1},
    {"bitarray_getBitStatus", (DL_FUNC)&bitarray_getBitStatus, 1},
    {"bitarray_Flip", (DL_FUNC)&bitarray_Flip, 1},
    {"createFile", (DL_FUNC)&createFile, 6},
    {"writeSlice", (DL_FUNC)&writeSlice, 4},
    {"readSlice", (DL_FUNC)&readSlice, 4},
    {NULL, NULL, 0}
};

void R_init_ncdfFlow(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);

}
