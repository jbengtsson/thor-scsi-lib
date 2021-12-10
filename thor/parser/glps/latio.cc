#include "latio.h"
#include "lattice.h"
#include "util.h"
#include "string.h"

#ifdef WITH_HDF5
#include "hdf5.h"
#endif

#include <vector>
#include <string>
#include <ctime>

#include <iostream>
#include <iomanip>

using namespace std;


size_t str_index(const vector<string>& vs, const string& s)
{
    for (size_t i = 0; i < vs.size(); ++i) {
        if (str_case_cmp(vs[i].c_str(), s.c_str()) == 0 )
            return i;
    }
    return -1;
}

size_t str_index(const vector<string>& vs, const char* s)
{
    for (size_t i = 0; i < vs.size(); ++i) {
        if (str_case_cmp(vs[i].c_str(), s) == 0 )
            return i;
    }
    return -1;
}


#ifdef WITH_HDF5
int write_rec(hid_t latgrp, const SP_STMT_LST* stmt)
{
    double dbasic[] = {0.0, 0.0};
    const hsize_t RANK = 1;
    hsize_t ddim[] = {1};
    ddim[0] = sizeof(dbasic)/sizeof(double);
    // cout << "dim: " << ddim[0] << ' ' << stmt->mag <<'|' << strlen(stmt->mag) << endl;

    //hid_t atype;
    hid_t dtype, dspace, dset, attr;
    herr_t ret;

    ddim[0] = 1;
    dspace = H5Screate(H5S_SCALAR); //H5Screate_simple(RANK, ddim, NULL);
    dtype  = H5Tcopy(H5T_C_S1);
    H5Tset_size(dtype, strlen(stmt->mag));
    dset = H5Dcreate(latgrp, stmt->name, dtype, dspace, H5P_DEFAULT);
    ret = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, (stmt->mag));

    SP_PRPT_LST* prpt = stmt->property;
    double *d;
    while(prpt) {
            if (prpt->nval > 0) {
                d = new double[prpt->nval];
                for (size_t i = 0; i < prpt->nval; ++i) {
                    d[i] = prpt->val[i];
                }
                ddim[0] = prpt->nval;
                ret = H5Sset_extent_simple(dspace, RANK, ddim, NULL);
                attr = H5Acreate(dset, prpt->property, H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT);
                ret = H5Awrite(attr, H5T_NATIVE_DOUBLE, d);
            } else {
                dspace = H5Screate(H5S_SCALAR);
                H5Tset_size(dtype, strlen(prpt->str));
                attr = H5Acreate(dset, prpt->property, dtype, dspace, H5P_DEFAULT);
                ret = H5Awrite(attr, dtype, prpt->str);
            }
        
        prpt = prpt->next;
    }
    /*    hid_t datatype = 
    hid_t elem = H5Dcreate(latgrp, "element", datatype, dataspace,
			H5P_DEFAULT);
    */
    return 0;
}

void write_lat_h5(const char* filename, const char* groupname, const SP_STMT_LST* stmt)
{
    hid_t       file, dataset;         /* file and dataset handles */

    FILE* pf = fopen(filename, "r");
    if (!pf) {
        file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }else {
        file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    
    char buf[128];
    string grp(groupname);
    time_t t = time(NULL);
    struct tm* local = localtime(&t);
    strftime(buf, 127, "_%Y%m%d_%H%M%S", local);
    grp += "_" + random_hexstring(4);
    grp.append(buf);

    //dataset = H5Gopen(file, groupname);
    //if (dataset >= 0) {
    //    H5Gunlink(file, groupname);
    //} 

    dataset = H5Gcreate(file, grp.c_str(), 0);

    const SP_STMT_LST* pstmt = stmt;
    while(pstmt) {
        write_rec(dataset, pstmt);
        pstmt = pstmt->next;
    }

    H5Gclose(dataset);
    H5Fclose(file);    
}
#endif


#ifdef DEBUG
int main(int argc, char* argv[])
{
    parse_lattice("../test/test.lat");

    STMT_LST *pstmt = get_statement_list();
    
    ios_base::fmtflags stdfmt = cout.flags();

    cout << endl
         << "## ------------------------------------------------------- ##" << endl
         << "## Lattice: test.lat "<< endl
         << "## ------------------------------------------------------- ##" << endl
         << endl;

    AP_LAT_REC rec[1000];
    int irec = 0;

    vector<string> element, magtype, property, vstring;
    int ielement, imagtype, iproperty, ivstring;

    hid_t       file, dataset;         /* file and dataset handles */

    file = H5Fcreate("latio.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dataset = H5Gcreate(file, "/lattice", 0);
    
    while(pstmt) {
        // loop over all statements
        cout << "[ " << pstmt->name << " : " << pstmt->mag
             << "  " << pstmt->type << "  " << pstmt->ifam
             << " ]" << endl;

        if((ielement = str_index(element, string(pstmt->name))) < 0) {
            element.push_back(string(pstmt->name));
            ielement = element.size() - 1;
        }

        if ((imagtype = str_index(magtype, pstmt->mag)) < 0) {
            magtype.push_back(string(pstmt->mag));
            imagtype = magtype.size() - 1;
        }

        /* */
        write_rec(dataset, pstmt);

        PRPT_LST *prpt = pstmt->property;

        while (prpt) {

            if((iproperty = str_index(property, prpt->property)) < 0) {
                property.push_back(string(prpt->property));
                iproperty = property.size() - 1;
            }


            cout << "    " 
                 << prpt->property << " =";
            
            if (prpt->nval > 0) {
                for (size_t i = 0; i < prpt->nval; ++i) {
                    cout << ' ' << setprecision(10) << prpt->val[i];
                    //
                    rec[irec].element = ielement;
                    rec[irec].magtype = imagtype;
                    rec[irec].property = iproperty;
                    rec[irec].vlength = prpt->nval;
                    rec[irec].vindex  = static_cast<int>(i);
                    rec[irec].vdouble = prpt->val[i];

                    rec[irec].vint = rec[irec].vstring = -1;
                    ++irec;
                }
                cout << endl;
            } else {
                // a string type
                cout << ' ' << prpt->str << endl;
                if((ivstring = str_index(vstring, prpt->str)) < 0) {
                    vstring.push_back(string(prpt->str));
                    ivstring = vstring.size() - 1;
                }

                //
                rec[irec].element = ielement;
                rec[irec].magtype = imagtype;
                rec[irec].property = iproperty;
                rec[irec].vlength = 0;
                rec[irec].vindex  = 0;
                rec[irec].vdouble = 0;
                rec[irec].vint = 0;
                rec[irec].vstring = ivstring;
                ++irec;
            }
            cout.flags(stdfmt);
            prpt = prpt->next;
        }
        cout << endl;
        pstmt = pstmt->next;
    }
    int N = irec;
    free_lattice();

    
    H5Gclose(dataset);
    H5Fclose(file);

    return 0;
}

#endif
