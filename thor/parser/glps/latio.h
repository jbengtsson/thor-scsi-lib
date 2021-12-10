#include "lattice.h"

#ifndef AP_LATIO_H_
#define AP_LATIO_H_

#ifdef __cplusplus
extern "C" {
 #endif 

/*
// lattice record.
struct AP_LAT_REC {
    int element;        // index of element name
    int magtype;        // index of magnet type
    int property;       // index of property
    int line;           // line number in original lattice
    int vtype;          // index of value type
    int vindex;         // index if this is part of an array
    int vlength;        // length if this is part of an array
    int vstring;        // index of string
    int vint;           // value
    double vdouble;     // value
};
*/

void write_lat_h5(const char* filename, const char* groupname, const SP_STMT_LST* stmt); 

#ifdef __cplusplus
};
 #endif 

#endif
