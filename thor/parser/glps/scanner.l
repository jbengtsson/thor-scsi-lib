/* -*- Mode: c; tab-width:4; indent-tabs-mode: nil; c-basic-offset:4 -*-  */
%{
/* $Id: scanner.l 386 2009-09-07 22:19:03Z doudoubai $
 *
 * Copyright (C) 2008 Lingyun Yang.
 *
 * This software may be used and distributed according to the terms
 * of the GNU General Public License, incorporated herein by reference.
 *
 * For more information please contact lyyang@lbl.gov
 *
 */


/* Author: L. Yang, lyyang@lbl.gov */
/* $LastChangedRevision: 386 $ */
/* $LastChangedDate: 2009-09-07 18:19:03 -0400 (Mon, 07 Sep 2009) $ */

#include "lattice.h"
#include "parser.h"
#include <math.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE_BUF_SIZE 1024

char* element_name;
unsigned int lineno = 1;

char linebuf[MAX_LINE_BUF_SIZE];
char strbuf[MAX_LINE_BUF_SIZE];

extern char** beam_lines;
extern int yydebug;
extern SP_SYMB_LST* symtab;
/* open multiple lattice file, push/pop the stack. */
#define MAX_INCLUDE_DEPTH 10
static int lat_file_stack_ptr = 0;
static YY_BUFFER_STATE lat_file_stack[MAX_INCLUDE_DEPTH]; 
static int lat_file_lineno[MAX_INCLUDE_DEPTH];

 /* processing property statement in an definition 
  * e.g. L = 0.1, property=(1,2,3,4,5); 
  * every time, when meet an identifier, we record its name,
  * in case we need it on the right hand side. we have to keep its value later.
  * i.e. vec/vec_size/value 
  */
#define NEW_PROPERTY_SYMBOL(stxt)   \
 yylval.symp = (SP_SYMB_LST*) malloc(sizeof(SP_SYMB_LST)); \
 yylval.symp->funcptr = NULL; \
 yylval.symp->name = strdup(stxt);                 \
 yylval.symp->vec_size = 0;                     \
 yylval.symp->vec = NULL;                       \
 yylval.symp->s = NULL;                       \
 yylval.symp->value = 0.0; \
 yylval.symp->next = NULL;


%}

  /* definitions */

  /* NUM   -?[0-9]+ */
  /* REAL  -?(([0-9]+)|(([0-9]*\.)?[0-9]+)([eE][-+]?[0-9]+)?)  */
NUM   [0-9]+
REAL  (([0-9]+)|(([0-9]*\.)?[0-9]+)([eE][-+]?[0-9]+)?) 
WS    [ \t\n]*
IDENTIFIER  [A-Za-z][A-Za-z0-9_]*

  /* x exclusive start-condition: */
%x sc_comment sc_str sc_beamline sc_element sc_action

%%

  /* --------------------------------------- */

  /* rules */
<INITIAL><<EOF>>  { 
 #ifdef DEBUG
    fprintf(stderr, "--- Reach the End of Lattice without an \"END\" statement. --\n"); 
 #endif

    if ( --lat_file_stack_ptr < 0) {
        yyterminate();
    }else {
        fprintf(stderr, "Switching...\n");
        yy_delete_buffer(YY_CURRENT_BUFFER);
        yy_switch_to_buffer(lat_file_stack[lat_file_stack_ptr]);
        lineno = lat_file_lineno[lat_file_stack_ptr];
    }
}

<sc_element,sc_action,sc_beamline><<EOF>>  {
    fprintf(stderr, "Reach the end of lattice file,"
            "but the statement is not finished with ;\n");
    return 0;
     }

  /* Syntax error at first line will also be copied into linebuf */
<*>\n.*  { 
    ++lineno; 
    strncpy(linebuf, yytext+1, MAX_LINE_BUF_SIZE-1);
    linebuf[MAX_LINE_BUF_SIZE-1] = '\0';
    yyless(1); 
 #ifdef DEBUG
    fprintf(stdout, "%d: %s\n", lineno, linebuf); 
 #endif
   }

   /* Comment */ 
"{"                    { BEGIN(sc_comment); }
<sc_comment>[^"}"\n]+  ;
<sc_comment>\n         {++lineno;}
<sc_comment>"}"+       { 
    BEGIN(INITIAL); 
    /* fprintf(stderr, "--- COMMENT ---\n"); */}

!.*\n  { ++lineno; /* fprintf(stderr, "\n---- Comment ----\n"); */} 

<*>\"[^\"\n]*(\n|\")  { 
    /* String, doesn't span > 1 line */
    /* TODO: remove this constraint */
    /* assuming "string", here remove the " */
    yylval.str = (char* )malloc(strlen(yytext));
    strncpy(yylval.str, yytext+1, strlen(yytext)-2);
    yylval.str[strlen(yytext)-2] = '\0';
    return STRING; 
   }

<*>{REAL} { yylval.dval = atof(yytext); return NUMBER; }
<*>{NUM}+ { yylval.dval = atof(yytext); return NUMBER; }

include[ \t\n]*","[ \t\n]*[^;\t\n]+";" {
    char *head, *tail, *inc, *incfile;
    size_t sz;
    if (lat_file_stack_ptr >= MAX_INCLUDE_DEPTH) {
        fprintf(stderr, "Includes nested too deeply\n");
        exit(-1);
    }

    fprintf(stderr, "Include file: %s\n", yytext);
    head = strchr(yytext, '"');
    tail = strrchr(yytext, '"');
    sz = strspn(head, tail);
    inc = (char*) malloc(sz*sizeof(char));
    incfile = inc;
    ++head;
    while(head != tail){
        *inc = *head;
        ++inc;
        ++head;
    }
    *inc = '\0';
    /* fprintf(stderr, "File: %s\n", incfile); */

    lat_file_lineno[lat_file_stack_ptr] = lineno;
    lat_file_stack[lat_file_stack_ptr++] = YY_CURRENT_BUFFER;

    yyin = fopen(incfile, "r");

    if (!yyin) {
        fprintf(stderr, "Error: Cannot open file %s\n", incfile);
        exit(-1);
    }
    
    yy_switch_to_buffer(yy_create_buffer(yyin, YY_BUF_SIZE));
    lineno = 0;

    free(incfile);
   }

<INITIAL>end{WS}";" { /* logical end of file */
    #ifdef DEBUG
        fprintf(stderr, "-- End of Lattice --\n"); 
    #endif
    return 0; } 

<INITIAL>{IDENTIFIER}/{WS}":"{WS}"line"{WS}"=" {
    SP_SYMB_LST* sp;
    BEGIN(sc_beamline); 
    /* fprintf(stderr, "\n----> Line: %s\n", yytext); */
    /* NEW_PROPERTY_SYMBOL(yytext); */
    sp = sym_find_append(yytext);
    sp->lineno = lineno;
    yylval.symp = sp;
    
    return ID;}


<INITIAL>{IDENTIFIER}/{WS}: {
    BEGIN(sc_element);
    /* fprintf(stderr, "\n---> element: %s\n", yytext); */

    yylval.str = strdup(yytext);
    /* yylval.str = strdup("hello"); */

    return ELEMENT;
 }

<INITIAL>{IDENTIFIER}/{WS},  { 
    BEGIN(sc_action); 
    /* fprintf(stderr, "Action: %s\n", yytext); */
    yylval.str = strdup(yytext);
    return ACTION; 
 }

<sc_beamline>LINE/{WS}=   { 
    /* fprintf(stderr, "\n---> beamline: \n", yytext); */
    return LINE; }

<sc_action>{IDENTIFIER}/{WS}= {
    /* convert_key_case(yytext); */
 #ifdef DEBUG
    fprintf(stderr, "Action: %s\n", yytext);
 #endif  
    
    NEW_PROPERTY_SYMBOL(yytext);
    return MAG_PROPERTY; }

<sc_element>{IDENTIFIER}/{WS}= { 
    /* general element */ 
 #ifdef DEBUG
    fprintf(stderr, "Magnet Property: %s\n", yytext);
 #endif

    NEW_PROPERTY_SYMBOL(yytext);
    return MAG_PROPERTY; 
}

<sc_beamline>{IDENTIFIER} { 
    /* BL: Line=(ID1, ID2, ...) */
    int i = beamline_def(yytext);
    /* fprintf(stdout, "Define beam line: %s\n", yytext); */
    char* ret;
    if (i > 0) {

     #ifdef DEBUG
        fprintf(stdout, "[DBG] Define beam line: %s\n", yytext);
     #endif

        yylval.str = strdup(beam_lines[i]);
    }else {
        if (!statement_def(yytext)) {
            fprintf(stderr, "Line: %d Element: |%s| is not defined\n", lineno+1, yytext);
            exit(-1);
        }

     #ifdef DEBUG
        fprintf(stdout, "[DBG]        element: %s\n", yytext);
     #endif

        yylval.str = strdup(yytext);
    }
    return BL; 
 }

<*>{IDENTIFIER} { 
    SP_SYMB_LST* sp = sym_find_append(yytext);
    yylval.symp = sp;
  
 #ifdef DEBUG
    fprintf(stdout, "\n[DBG] Got ID: %s\n", yytext);  
 #endif
    
    return ID; 
}


<*>{IDENTIFIER}"."{IDENTIFIER} { 
    /* NEW_PROPERTY_SYMBOL(yytext); */
    SP_SYMB_LST *sp;
    SP_STMT_LST *pstmt;
    SP_PRPT_LST *prpt;
    char *p, *pelem, *prop, *s;
    double d;
    size_t i;
    p = strdup(yytext);
    pelem = p;
    while(*p != '.') ++p;
    *p = '\0';
    prop = p+1;

    s = NULL;
    
    /* searching for element and its property */
    pstmt = statement_def(pelem);
    if (!pstmt) {
        fprintf(stderr, "ERROR: element %s not defined\n", pelem);
        exit(-1);
    }

    prpt = get_statement_property(pstmt, prop);

    if (!prpt) {
        fprintf(stderr, "ERROR: element %s has no property %s defined\n", pelem, prop);
        exit(-1);
    }

    sp = sym_find_append(yytext);
    i = 0;
    /* if (sp->name) fprintf(stderr, "name : %s\n", sp->name); */
    if (prpt->str) sp->s = strdup(prpt->str);
    if (prpt->nval == 1) sp->value = prpt->val[0];
    else if (prpt->nval > 1) {
        sp->vec_size = prpt->nval;
        sp->vec = (double *)malloc(prpt->nval*sizeof(double));
        for (i = 0; i < sp->vec_size; ++i) {
            sp->vec[i] = prpt->val[i];
        }
    }
    yylval.symp = sp;
    free(pelem);

    return ID;
}


<sc_beamline>"-" { return INV; }

<*>"+"|"-"|"*"|"/" { return yytext[0]; }

<*>":"|","|"="|"("|")" { return yytext[0]; }

<*>";" { BEGIN(INITIAL); return yytext[0]; }

<*>[ \t\r] ;

<*>"&" ;

. {fprintf(stderr, "\nUnknown:|%c|, Skip\n", yytext[0]); }

%%

