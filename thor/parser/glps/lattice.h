/*
 * $Id: lattice.h 386 2009-09-07 22:19:03Z doudoubai $
 *
 * Copyright (C) 2008 Lingyun Yang.
 *
 * This software may be used and distributed according to the terms
 * of the GNU General Public License, incorporated herein by reference.
 *
 * For more information please contact lyyang@lbl.gov
 *
 */


/* Author: Lingyun Yang, lyyang@lbl.gov */

/* $LastChangedRevision: 386 $ */
/* $LastChangedDate: 2009-06-07 19:16:19 -0700 (Sun, 07 Jun 2009) $ */

#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef LATTICE_H
#define LATTICE_H

#ifdef __cplusplus
extern "C" {
 #endif 

    #define GLPS_VER_MAJOR 1
    #define GLPS_VER_MINOR 0
    #define GLPS_VER_PATCH 2
    #define GLPS_VER_STR "$Rev: 386 $"

    /*! This defines the type of elements, it must be same as Tracy, otherwise,
     *  the output flat machine file would be different.
     * T_MARKER = -1, , T_DRIFT=0, T_MULTIPOLE=1,
     * T_CAVITY= 2,  T_THINKICK=  3,  T_WIGGLER = 4}
     */
    enum SP_ELEMENT_TYPE {
        ELEMT_USERDEF = 0, ELEMT_MARKER    =-1, 
        ELEMT_DRIFT   = 2, ELEMT_MULTIPOLE =3,
        ELEMT_CAVITY  = 4, ELEMT_WIGGLER   =5,
        ELEMT_THINKICK= 6,
        STMT_LINE = 64,
        STMT_ACTION_SET = 128, STMT_ACTION_TITLE          
    } ;


    enum SP_VALUE_TYPE {T_DOUBLE = 0, T_DOUBLE_ARRAY,  T_STRING};

    /*! Global symbol table, store every variable, function defined. It also
     *  stores the property of magnets and actions. But since the magnet
     *  properties are copied into a new structure after matched in Yacc
     *  grammar, and freed after that, they are seperated from the global
     *  linked list.
     */
    struct symbol_table_rec {
        struct symbol_table_rec *next;
        double     (*funcptr)();    /* cos, sin, ... */
        char              *name;    /* PI, L, ... */
        /* char      *element;*/    /* PUE01, ..... */
        /* char       *family;*/    /* Quad, BPM, Bend, Sextupole, ... */
        /* unsigned int  ifam;*/    /* Family number, the order read in*/
        size_t           lineno;    /* line number of this symbol */
        size_t         vec_size;    /* */
        double             *vec;    /* */
        double            value;    /* */
        char                 *s;    /* */
    };
    typedef struct symbol_table_rec SP_SYMB_LST;

    /*! This stores the property of a magnet or action. It is used by other
     *  data structures.
     */
    struct property_table_rec {
        struct property_table_rec *next;
        char    *property;  /* */
        size_t       nval;  /* */
        double       *val;  /* */
        char         *str;
    };
    typedef struct property_table_rec SP_PRPT_LST;

    /*! A linked list for magnets. It uses property tables for each
     *  element. After a maching in YACC grammar, it gets data from a local
     *  symbol table, and release the table. When refering the property of a
     *  magnet, e.g. QF.L, We have to check both the global symbol table, and
     *  this magnet table, and keep them synchronized.
     */
    struct statement_rec{
        struct statement_rec *next;
        char* name;
        char* mag;
        int   type;
        int   ifam;       /* Family Number, the order defined */
        SP_PRPT_LST *property;
    };
    typedef struct statement_rec SP_STMT_LST;

    /*!\brief parse the lattice file into internal data structure 
     */
    int parse_lattice(const char* f);

    /*! \brief release the internal structure storing parsed lattice 
     */
    int free_lattice();

    char *get_path_root(char *f);
    int parse_lattice_flat(char* f, char* flat);
    int parse_lattice_xml(char* lat, char* xml);
    int parse_lattice_ini(char* f);
    int print_flat_lattice(FILE* pf, char* bl);


    /* ----------------------------------------- */

    SP_SYMB_LST *sym_find_append(char* s);
    SP_SYMB_LST* sym_remove_id(char* s);    
    SP_SYMB_LST* sym_remove_node(SP_SYMB_LST* s);    

    /*! \brief get the definition of a statement, including action, element
     *   and beamline 
     */
    SP_STMT_LST *statement_def(const char* elem);
    int beamline_def(char* bl);
    int save_elements_db(SP_STMT_LST* pelem, char* dbfile);
    int save_lat_db(SP_STMT_LST* pelem, const char* dbfile);
    /*! \brief compare two string, case-incensitive */
    int str_case_cmp(const char* ps1, const char* ps2);

    /*! \brief Print symbol table */
    void print_symb_tab(SP_SYMB_LST* p);
    /*! \brief Print statement from pelem */
    int print_elements(SP_STMT_LST* pelem);

    SP_STMT_LST* get_statement_list();
    SP_PRPT_LST *get_statement_property( SP_STMT_LST* elem, char* property);
    int get_statement_property_double(const char* statement, const char* property, double* val);
    int get_statement_property_text(const char* statement, const char* property, char** val);

    /*
    struct keyword_dic{
        int id;
        char* class;
        char* property;
    } KEYWORD[] = {
        {0, "MARKER,MARK", ""}, 
        {1, "DRIFT,DRIF", "L"}};

    enum {MARKER=0, BPM = 1,
          DRIFT=2, DRIF=2,
          QUAD=3, QUADRUPOLE = 3,
          SEXTUPOLE=4};
    */

    void get_parser_version(int* vmaj, int* vmin, int* vpat, int* vrev);

#ifdef __cplusplus
}
#endif

#endif
