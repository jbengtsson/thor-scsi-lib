/* -*- Mode: c; tab-width:4; indent-tabs-mode: nil; c-basic-offset:4 -*-  */
%{
 /*
  * $Id: parser.y 386 2009-09-07 22:19:03Z doudoubai $
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
 /* $LastChangedDate: 2009-09-07 18:19:03 -0400 (Mon, 07 Sep 2009) $ */

 #include "lattice.h"
 #include "latio.h"
 #include <stdio.h>
 #include <string.h>
 #include <math.h>
 #include <ctype.h>

 #ifndef M_PI
 #define M_PI 3.14159265358979323846
 #endif

 #define XML 1

 /* int yydebug=1; */

 /*! global symbol table (linked list) */
 /*  static SP_SYMB_LST* symtab = NULL; */
 SP_SYMB_LST* symtab = NULL;

 extern unsigned int lineno;
 extern char linebuf[];
 extern char strbuf[];

 /*! lattice title */
 static char* lat_title = NULL;

 /*! element table (linked list) */
 static SP_STMT_LST *stmt_table = NULL;

 /*! tail to element table (linked list), speed up appending new element */
 static SP_STMT_LST *stmt_table_tail = NULL;

 static int echo = 1;
 /*! Flat file */
 /* static FILE* pflat; */

 /* ["key_0", "elements_0", "key_1", "elements_1", ... ] */
 char** beam_lines; 

 /*! size of beam_lines array, twice of the number of beam lines */
 int n_beam_lines = 0;


 /* if yywrap returns 0, the scanner continues scanning, while if it returns 1,
   the scanner returns a zero token to report the end of file. If yywrap()
   returns 0 to indicate that there is more input, it needs first to adjust
   yyin to point to a new file, probably using fopen()
 */
 extern FILE* yyin;
 int yywrap(void) { 
     return 1; 
 };
 int yylex();
 int yyerror(char* s);
 

 /*! Free symbol table (linked list) */
 static void free_symb_tab(SP_SYMB_LST* head);

 void sym_addfunc(char* name, double(*func)(double));
 SP_SYMB_LST *sym_append();
 
 /* SP_STMT_LST* statement_def(char* elem); */
 SP_STMT_LST* statement_add(char* elem);
 void statement_add_property_list(SP_STMT_LST* elem, 
                                SP_SYMB_LST *head);
 
 char* beamline_dup(char* bl);
 void beamline_add(char* bl, char* elem);
 char* beamline_reverse(char* bl);
  
 void free_stmt_table();
 void free_beamline();
  
%}

%union {
    double  dval;   /* plain number */
    char*    str;   /* string */
    SP_SYMB_LST *symp;     /* variable */
    struct symbol_pt_list *symp_list;
}

%token <str> STRING ELEMENT ACTION
%token <dval> REAL
%token <dval> NUMBER
%token <symp> MAG_PROPERTY
%token <symp> MAG_PROPERTY_VEC
%token <symp> ID ID_PROPERTY
%token <dval> FUNC

%token <str> BL

%token SET TITLE SHOW INCLUDE
%token BEND BPM CAVITY CORR DRIFT MARKER MULTIPOLE INV LINE QUAD SEXT WIGGLER 

%left ','
%left '-' '+'
%left '*' '/'
%nonassoc UMINUS

%type <dval> expression
%type <symp> expression_list
%type <symp> mag_property
%type <symp> mag_property_list
%type <str> beamline

%%

statement_list: statement 
        |       statement_list statement
        ;

statement: ID '=' expression ';' { 
    $1->value = $3;
    /* printf("%s = %f\n", $1->name, $3); */
   }
|  ID '=' '(' expression_list ')' ';' { 
    SP_SYMB_LST *p = $4, *q;
    size_t i=0, n = 0;
    /* printf(" %s:\n", $1->name); */
    while(p) {
        ++n;
        /*
          if (p->s) printf("    %s\n", p->s);
          printf("    %f %d (", p->value, p->vec_size);
          for (i = 0; i < p->vec_size; ++i) printf(" %f", p->vec[i]);
          printf(")\n");
        */
        p = p->next;
    }
    $1->vec_size = n;
    $1->vec = (double*) malloc(n*sizeof(double));
    p = $4;
    while (p) {
        $1->vec[i++] = p->value;
        q = p->next;
        free(p->s);
        free(p);
        p = q;
    }
   }
|ACTION ',' ID ';' {
    size_t i =0;
    if (str_case_cmp($1,"PRINT") == 0) {
        if($3->name) fprintf(stdout, "%s : ", $3->name);
        if($3->vec) {
            for(i=0; i< $3->vec_size; ++i) {
                fprintf(stdout, " %f", $3->vec[i]);
            }
        }
        if ($3->s) fprintf(stdout, "%s", $3->s);
        if($3->value) fprintf(stdout, "%f", $3->value);
        
        fprintf(stdout, "\n");
    }
    free($1);
 }
|ACTION ',' STRING ';' { 
    SP_STMT_LST *p = statement_add($1);
    /* strncpy(lat_title, strbuf+1, strlen(strbuf)-2);*/
    if(lat_title) free(lat_title); 
    /* lat_title = strdup($3);  */
    lat_title = (char*)malloc(strlen($3)+1);
    strcpy(lat_title, $3);

    p->type = STMT_ACTION_TITLE;
    if (p->mag) free(p->mag);
    /* p->mag = strdup($1); */
    p->mag = (char*) malloc(strlen($1)+1);
    strcpy(p->mag, $1);

    p->ifam -= 1; /* do not increase the family index, since this is an action */
    p->property = (SP_PRPT_LST*) malloc(sizeof(SP_PRPT_LST));
    p->property->next = NULL;
    p->property->property = strdup("TITLE");
    p->property->nval = 0;
    p->property->val = NULL;
    p->property->str = strdup(lat_title);

    if (echo) {
        fprintf(stdout, "[ %s : %s ]\n", $1, $3);
    }
    free($1);
    free($3);
   }
|  ACTION ',' mag_property_list ';' {
    SP_STMT_LST *p =  statement_add($1); 
    /* fprintf(stdout, "--> mp %s %f %s\n", $3->name, $3->value, $3->s); */
    /* fprintf(stdout, "--> mp %s %f %s\n", $3->next->name, $3->next->value, $3->next->s); */
    statement_add_property_list(p, $3);
    p->type = STMT_ACTION_SET;
    p->mag = strdup($1);
    p->ifam -= 1; /* do not increase the family index, since this is an action */
    /* act_setup_beam($3); */
    /* printf("-------------------%s-------------\n", $1);
       print_symb_tab(symtab); */
    /* messed up str in symtab if uncomment the following two lines */
    if (str_case_cmp($1, "save") == 0) {
        char *filename, *groupname;
        get_statement_property_text("set", "output", &filename);
        get_statement_property_text("save", "name", &groupname);
        #ifdef WITH_HDF5
        write_lat_h5(filename, groupname, stmt_table);
        #endif
    }
    free($1);
    free_symb_tab($3);
   }
| ELEMENT ':' ID ';' { 
    SP_STMT_LST *p =  statement_add($1);    
    p->mag = strdup($3->name);
    p->type = ELEMT_THINKICK;
    /* the ID has been added to the symbtab in scanner, remove it here */
    sym_remove_id($3->name);
    /* printf("BPM: %s\n", $1->name); */}
| ELEMENT ':' ID ',' mag_property_list ';' {
    SP_STMT_LST *p =  statement_add($1); 
    statement_add_property_list(p, $5);
    p->type = ELEMT_CAVITY;
    p->mag = strdup($3->name);
    free_symb_tab($5);
    /* free_symb_tab($3); */
    free($1);
    /* printf("Element: %s %s\n", $1, $3->name); */
    /* the ID has been added to the symbtab in scanner, remove it here */
    sym_remove_id($3->name);
   }

|  ID ':' LINE '=' '(' beamline ')' ';'  {
    SP_STMT_LST *ps;
    SP_PRPT_LST *pp;
    /* printf("LINE: %s = ( %s )\n", $1->name, $6); */
    beamline_add($1->name, $6);
    ps = statement_add($1->name);
    ps->mag = strdup("BEAMLINE");
    ps->type = STMT_LINE;
    pp = get_statement_property(ps, "LINE");
    pp->str = strdup($6);
    /* sym_remove_node($1); */
    /* free_symb_tab($1); */
    /* p = symtab; */
    /* while(p != NULL && p != $1) { */
        /* fprintf(stderr, "%s %x\n", p->name, p); */
    /*    p = p->next;
          }*/
    /* if (echo && p == $1) {fprintf(stdout, "Found\n");}
       else if (echo && p == NULL) {fprintf(stderr, "Not found\n");} */
    $1->s = strdup($6);
    free($6);
   }
|  ID ID ';' {
    /* such as old tracy lattice, "define lattice;" */
    fprintf(stderr, "ERROR: Line %d, Syntax error:\n  %s %s;\n", lineno, $1->name, $2->name);
    if (str_case_cmp($1->name, "define") == 0 && 
        str_case_cmp($2->name, "lattice") == 0) {
        fprintf(stderr, "\nThe old Tracy-II lattice is obsolete, please check the document for new format.\n");
    }
    exit(-1);
 }
;

beamline: BL { 
    $$ = $1;
    /* fprintf(stderr, "---+-> %s\n", $$); */
 }
|  beamline ',' beamline {
    char* ret = (char*) malloc(sizeof(char)*(strlen($1) + strlen($3) + 2));
    strcpy(ret, $1);
    strcat(ret, ",");
    strcat(ret, $3);
    /* fprintf(stderr, "---+-> %s + %s = %s\n", $1, $3, ret); */
    free($1);
    free($3);
    $$ = ret;
   }
|  NUMBER '*' '(' beamline ')' {
    int i = 0;
    int n = strlen($4)+1;
    int ndup = (int)($1);
    char *ret, *stmp;
    ndup = abs(ndup);

    ret = (char*) malloc(sizeof(char)*n*(ndup+1));
    stmp = ($1<0.0) ? beamline_reverse($4) : strdup($4);

    strcpy(ret, stmp);
    for (i = 1; i < ndup; ++i) {
        strcat(ret, ",");
        strcat(ret, stmp);
    }
    free($4);
    free(stmp);
    $$ = ret;
   }
|  NUMBER '*' BL {
    int n = strlen($3)+1; /* need 1 more space for "," */
    int ndup = $1, i = 0;
    char *ret, *stmp;
    ndup = abs(ndup);
    ret = (char*) malloc(sizeof(char)*n*(ndup+1));

    if ($1 < 0) stmp = beamline_reverse($3);
    else stmp = strdup($3);

    /* fprintf(stderr, "---+-> %d * %s, n = %d\n", ndup, $3, n*(ndup+1)); */
    /* fprintf(stderr, "---+-> %d * %s, n = %d\n", ndup, stmp, n*(ndup+1)); */
    strcpy(ret, stmp);
    
    for (i = 1; i < ndup; ++i) {
        /* fprintf(stderr, "%d: %s  %s\n", i, ret, stmp); */
        strcat(ret, ",");
        strcat(ret, stmp);
    }

    free($3);
    free(stmp);
    /* fprintf(stderr, "---+-> %s\n", ret); */
    $$ = ret;
   }
|  INV BL { 
    char* ret;
    /* fprintf(stderr, "Reverse Beam Line: %s\n", $2); */
    ret = beamline_reverse($2);
    free($2);
    $$ = ret;
 }
|  INV '(' beamline ')' {
    char* ret;
    /* fprintf(stderr, "Reverse Beam Line: %s\n", $3); */
    ret = beamline_reverse($3);
    free($3);
    $$ = ret;
   }
;

mag_property_list: mag_property { 
    $$ = $1;
 }
| mag_property_list ',' mag_property {
    SP_SYMB_LST *p = $1;
    $$ = $1;
    while(p->next) p = p->next;
    p->next = $3;
}
;
  
mag_property: MAG_PROPERTY '=' expression {
    /* */
    $1->value = $3;
    if ($1->next) {
        fprintf(stderr, 
                "--+-> PROPERTY should be independent!\n");
        fprintf(stderr, "  %s = %f\n", $1->name, $3);
        exit(-1);
    }
    $$ = $1;
    $$->s = NULL;
    $$->next = NULL;
}
| MAG_PROPERTY '=' STRING {
    $1->s = strdup($3);
    free($3);
    if ($1->next) {
        fprintf(stderr, 
                "--+-> PROPERTY should be independent!\n");
        fprintf(stderr, "  %s = %s\n", $1->name, $3);
        exit(-1);
    }
    $$ = $1;
    $$->vec_size = 0;
    $$->vec = NULL;
    $$->value = 0;
    $$->next = NULL;
    /* fprintf(stdout, "STRING 1: %s = %s %d\n", $$->name, $$->s, strlen($$->s)); */
    /* fprintf(stdout, "STRING 1: %s = %s %d\n", $1->name, $3, strlen($3)); */
  }
| MAG_PROPERTY '=' '(' expression_list ')' {
    SP_SYMB_LST *p = $4, *q;
    /* printf("PropertyVec: %s = (", $1->name); */
    /* printf("%s = (", $1->name); */
    size_t i=0, n = 0;
    while(p) {
        ++n;
        /* printf("    %s %s %s = %f \n", p->element, p->family, p->name, p->value); */
        p = p->next;
    }
    $1->vec_size = n;
    $1->vec = (double*) malloc(n*sizeof(double));
    p = $4;
    while (p) {
        $1->vec[i++] = p->value;
        q = p->next;
        free(p);
        p = q;
    }

    $$ = $1;
    $$->s = NULL;
    $$->next = NULL;
    /* printf(")\n"); */
  }       
;

expression_list: expression ',' expression {
    $$=(SP_SYMB_LST*) malloc(sizeof(SP_SYMB_LST));
    $$->value = $1;
    $$->s = NULL;
    $$->vec_size = 0;
    $$->vec = NULL;
    $$->next = (SP_SYMB_LST*) malloc(sizeof(SP_SYMB_LST)); 
    $$->next->value = $3;
    $$->next->s = NULL;
    $$->next->vec_size = 0;
    $$->next->vec = NULL;
    $$->next->next = NULL;
}
| expression_list ',' expression {
    SP_SYMB_LST* p = $1;
    while(p->next) p = p->next;
    p->next = (SP_SYMB_LST*) malloc(sizeof(SP_SYMB_LST));
    p->next->value = $3;
    p->next->s = NULL;
    p->next->vec_size = 0;
    p->next->vec = NULL;
    p->next->next = NULL;
    $$=$1;
}
;

expression: expression '+' expression { $$ = $1 + $3; }
|  expression '-' expression { $$ = $1 - $3; }
|  expression '*' expression { $$ = $1 * $3; }
|  expression '/' expression {
        if($3 == 0.0){
            yyerror("Divide by zero");  
        /* printf("%d\n", $3); */
        }
        else $$ = $1 / $3;
    }
|  '-' expression %prec UMINUS { $$ = -$2; }
|  '(' expression ')' { $$ = $2; }
|  NUMBER { $$ = $1; }
|  REAL { $$ = $1; }
|  ID { $$ = $1->value; }
|  ID_PROPERTY { $$ = $1->value; }
|  ID '(' expression ')' { /* for function like sin(pi) */
    if ( $1->funcptr ) {
      /* fprintf(stderr, "%s, %f %f\n", $1->name, $3, ($1->funcptr)($3)); */
      $$ = ($1->funcptr)( $3 );
    }
        else {
            printf("%s not a function\n", $1->name);
            $$ = 0.0;
        }
    }
;

%%

/** \brief Get version of parser library
 */
void get_parser_version(int* vmaj, int* vmin, int* vpat, int* vrev)
{
    char *s = strdup(GLPS_VER_STR);
    char *p = s + strlen(s) - 1;

    if (vmaj == NULL || vmin == NULL || vpat == NULL || vrev == NULL)
        return;

    *vmaj = GLPS_VER_MAJOR;
    *vmin = GLPS_VER_MINOR;
    *vpat = GLPS_VER_PATCH;
    *vrev = 0;

    while (!isdigit(*p) && p>= s) {
        *p = '\0';
        --p;
    }
    while(isdigit(*p) && p >= s) --p;
    if (p < s) return;
    
    *vrev = atoi(p+1);
    
    free(s);
}


/*!
 * \brief case-insensitive compare, similar return convention as strcmp.
 */
int str_case_cmp(const char* ps1, const char* ps2)
{
    char c1, c2;
    int v;
    do {
        c1 = *ps1++;
        c2 = *ps2++;
        if (c1 >='A' && c1 <= 'Z') c1 += 'a' - 'A';
        if (c2 >='A' && c2 <= 'Z') c2 += 'a' - 'A';        
        v = (int) c1 - (int) c2;
    }while ((v==0) && (c1 != '\0'));

    return v;
}


/*!
 *\brief check if symbol in the global linked list
 */
int get_sym_text(char* s, char** val)
{
    SP_SYMB_LST *sp = symtab;
    if (!symtab) {
        *val = NULL;
        return 0;
    }

    while(sp) {
        if (str_case_cmp(sp->name, s) == 0) {
            *val = strdup(sp->s);
            return 1;
        }
        if (!sp->next) break;
        else sp = sp->next;
    }

    return 0;
}


/*!
 *\brief check if symbol in the global linked list
 */
int get_sym_double(char* s, double* val)
{
    SP_SYMB_LST *sp = symtab;
    if (!symtab) {
        *val = 0.0;
        return 0;
    }

    while(sp) {
        if (str_case_cmp(sp->name, s) == 0) {
            *val = sp->value;
            return 1;
        }
        if (!sp->next) break;
        else sp = sp->next;
    }

    return 0;
}


/*!
 *\brief check if symbol in the global linked list
 */
int get_sym_double_vec(char* s, int* n, double** val)
{
    SP_SYMB_LST *sp = symtab;
    unsigned int i = 0;
    if (!symtab) {
        *val = NULL;
        *n = 0;
        return 0;
    }
    while(sp) {
        if (str_case_cmp(sp->name, s) == 0) {
            *val = (double*) malloc(sp->vec_size*sizeof(double));
            for (i = 0; i < sp->vec_size; ++i) (*val)[i] = sp->vec[i];
            *n = sp->vec_size;
            return 1;
        }
        if (!sp->next) break;
        else sp = sp->next;
    }

    return 0;
}

SP_SYMB_LST* sym_remove_id(char* s)
{
    /* char *p; */
    SP_SYMB_LST *head = symtab, *tail = symtab;
    if (!head) {
        return NULL;
    }

    if (str_case_cmp(symtab->name, s) == 0) {
        symtab = symtab->next;
        free(head->name);
        free(head->vec);
        free(head->s);
        free(head);
        return symtab;
    }

    head = symtab->next;
    while (head) {
        /* case insensitive */
        if (str_case_cmp(head->name, s) == 0) {
            tail->next = head->next;
            free(head->name);
            free(head->vec);
            free(head->s);
            free(head);
            return tail->next;
        }
        tail = head;
        head = head->next;
    }
    return NULL;
}


SP_SYMB_LST* sym_remove_node(SP_SYMB_LST* sp)
{
    /* char *p; */
    SP_SYMB_LST *head = symtab, *tail = symtab;
    if (!sp) {
        return NULL;
    }

    if (sp == symtab) {
        symtab = symtab->next;
        free(sp->name);
        free(sp->vec);
        free(sp->s);
        free(sp);
        return symtab;
    }

    head = symtab->next;
    while (head != sp && head) {
        tail = head;
        head = head->next;
    }

    if (head == sp) {
        tail->next = head->next;
        free(head->name);
        free(head->vec);
        free(head->s);
        free(head);
        return tail->next;
    }

    return NULL;
}

/*! 
 * \brief look up a symbol table entry, add if not present 
 *
 * \param s symbol, case insensitive.
 * \return An symbol table node.
 */
SP_SYMB_LST * sym_find_append(char* s)
{
    /* char *p; */
    SP_SYMB_LST *sp = symtab;
    if (!symtab) {
        symtab = (SP_SYMB_LST*) malloc(sizeof(SP_SYMB_LST));
        sp = symtab;
    } else {
        while (sp) {
            /* case insensitive */
            if (str_case_cmp(sp->name, s) == 0) return sp;
            if (!sp->next) break;
            else sp = sp->next;
        }
      /* printf("\nNot found %s\n", s); */
        sp->next = (SP_SYMB_LST*) malloc(sizeof(SP_SYMB_LST));
        sp = sp->next;
    }

    sp->name = strdup(s);
    sp->funcptr = NULL;
    /* sp->family = NULL;
       sp->element = NULL; */
    sp->next = NULL;
    sp->vec_size = 0;
    sp->vec = NULL;
    sp->value = 0;
    sp->s = NULL;
    return sp;
} /* symlook */

/*!
 * \brief Print symbol table
 */
void print_symb_tab(SP_SYMB_LST* p)
{
    unsigned int i = 0;
    SP_SYMB_LST* q = p;

    fprintf(stdout, "Functions: \n");
    while(q) {
        if (q->funcptr) {
            fprintf(stdout, " %s", q->name);
        }
        q = q->next;
    }
    fprintf(stdout, "\n");

    while (p) {
        if (p->funcptr) {
            p = p->next;
            continue;
        }

        fprintf(stdout, "%s", p->name);
        if (p->s) fprintf(stdout, " = \"%s\"\n", p->s);
        else if (p->vec_size > 0) {
            fprintf(stdout, "= [%zd] : (", p->vec_size);
            for (i = 0; i < p->vec_size; ++i) {
                fprintf(stdout, " %f", p->vec[i]);
            }
            fprintf(stdout, " )\n");
        }else {
            fprintf(stdout, " = %f\n", p->value);
        }
        p = p->next;
    }
    /* fprintf(stdout, "--------------- END -----------\n"); */
}

/*! \brief Free symbol table (linked list) 
 */
void free_symb_tab(SP_SYMB_LST* head)
{    
  SP_SYMB_LST *sp = head, *ph = head;
    while(ph) {
        sp = ph;
        ph = ph->next;
        /*
          if (sp->vec_size == 0)
          printf("sym: %s = %f\n", sp->name, sp->value);
          else {
          printf("sym: %s = ( ", sp->name);
          for(i = 0; i < sp->vec_size-1; ++i)
          printf("%f ", sp->vec[i]);
          printf("%f )\n", sp->vec[sp->vec_size-1]);
          }
        */

        if (sp->vec_size > 0) free(sp->vec);
        /* printf("name: address: %p\n", (void*) sp->name); */
        if (sp->s)     free(sp->s);
        if (sp->name)    free(sp->name);
        /* if (sp->element) free(sp->element); */
        /* if (sp->family)  free(sp->family); */

        free(sp);
    }
}

/*! \brief Add function of single variable into expressions.
 *
 * \param Function name, such as "sin", "cos", "tan".
 * \param pointer to the function, can be sin, cos, tan, and also a user defined function.
 *
 * \brief add function support, sin, exp, cos, .....
 */
void sym_addfunc(char* name, double (*func) (double))
{
    SP_SYMB_LST *sp = sym_find_append(name);
    sp->funcptr = func;
}

/*! \brief Error report with line number.
 *
 * If the error (syntax error)  happens at the first line, the position (line number) is not given correctly.
 */
int yyerror(char* s)
{
    fprintf(stderr, "Line %u: %s\n  %s\n", lineno, s, linebuf);
    exit(-1);
}


int get_statement_property_text(
    const char* statement, 
    const char* property, 
    char** val)
{
    SP_STMT_LST *st = stmt_table;
    while (st) {
        /* if (st->name) fprintf(stderr, "ST: %s\n", st->name); */
        if (str_case_cmp(st->name, statement) == 0) {
            /* fprintf(stderr, "  Yes\n"); */
            SP_PRPT_LST *pt = st->property;
            while (pt) {
                if (str_case_cmp(pt->property, property) == 0) {
                    *val = strdup(pt->str);
                    return 1;
                }
                pt = pt->next;
            }
        }
        st = st->next;
    }
    return 0;
}

int get_statement_property_double(const char* statement, const char* property, double* val)
{
    SP_STMT_LST *st = stmt_table;
    while (st) {
        if (str_case_cmp(st->name, statement) == 0) {
            SP_PRPT_LST *pt = st->property;
            while (pt) {
                if (str_case_cmp(pt->property, property) == 0) {
                    if (pt->nval > 0) *val = pt->val[0];
                    else *val = 0.0;
                    return 1;
                }
                pt = pt->next;
            }
        }
        st = st->next;
    }
    return 0;
}

/*! \brief Get the property table (linked list) of an element. If it does not
 * exist create one.
 *
 * \param elem pointer to an element.
 * \param propty property name.
 * \brief Get the property of an element.
 */
SP_PRPT_LST *get_statement_property(
    SP_STMT_LST *elem,
    char* propty)
{
    /* assuming elem is already there */
    SP_PRPT_LST *p = elem->property, *t = NULL;
    
    if (!p) {
        /* Add a new on, initialize */
        elem->property = (SP_PRPT_LST*) malloc(sizeof(SP_PRPT_LST));
        p = elem->property;
        p->str = NULL;
    } else {
        while (p != NULL && 
               str_case_cmp(propty, p->property) != 0) {
            t = p;
            p = p->next;
        }
        if (p == NULL) {
            /* Append a new one, initialize */
            t->next = malloc(sizeof(SP_PRPT_LST));
            p = t->next;
            p->str = NULL;
        }else{
            return p;
        }
    }
    
    p->property = strdup(propty);
    p->nval = 0;
    p->val = NULL;
    p->next = NULL;

    return p;
}

/*! \brief Add a property list to a element.
 *
 * append to the end.
 */
void statement_add_property_list(
    SP_STMT_LST *elem, 
    SP_SYMB_LST *head)
{       
    SP_SYMB_LST *p = head;
    int i = 0;
    
    SP_PRPT_LST *prpt;

    /* */
    while(p) {
        prpt = get_statement_property(elem, p->name);
        prpt->nval = 0;
        if (p->s) {
            prpt->str = strdup(p->s); 
            /* fprintf(stdout, "--> %s %s\n", prpt->property, prpt->str); */
        }else {
            prpt->str = NULL;
            if (prpt->val) free(prpt->val);
            if (p->vec_size > 0) {
                prpt->val = malloc(sizeof(double)*p->vec_size);
                for (i = 0; i < p->vec_size; ++i)
                    prpt->val[i] = p->vec[i];
                prpt->nval = p->vec_size;
            } else {
                prpt->val = malloc(sizeof(double));
                prpt->nval = 1;
                prpt->val[0] = p->value;
            }
        }
        p = p->next;
    }
}

/*! \brief Free the whole element table (linked list)
 */
void free_stmt_table()
{
    SP_STMT_LST *p = stmt_table;
    SP_PRPT_LST *prpt_p, *prpt_head;

    while(stmt_table) {
        p = stmt_table;
        stmt_table = stmt_table->next;
        free(p->name);
        if (p->mag) free(p->mag);
        prpt_head = p->property;
        while(prpt_head) {
            /* fprintf(stdout, "  %s %s\n", prpt_head->property, prpt_head->str); */
            prpt_p = prpt_head;
            prpt_head = prpt_head->next;
            
            free(prpt_p->property);
            if (prpt_p->val) free(prpt_p->val);
            if (prpt_p->str) free(prpt_p->str);
            free(prpt_p);
        }
        free(p);
    }    
}

/*!
 */
SP_STMT_LST* statement_def(const char* elem)
{
  SP_STMT_LST *p = stmt_table;
  
  while(p) {
    if (str_case_cmp(elem, p->name) == 0) return p;
    /* fprintf(stderr, "%s %s\n", elem, p->name); */
    p = p->next;
  }
  return NULL;
}

/*! \brief Add an element to the element table (linked list)
 *
 */
SP_STMT_LST* statement_add(char* elem)
{
    SP_STMT_LST *p = stmt_table;

    if (stmt_table == NULL) {
        stmt_table = malloc(sizeof(SP_STMT_LST));
        p = stmt_table;
        p->ifam = 1;
    }else {
        do {
            if (str_case_cmp(elem, p->name) == 0) return p;
            p = p->next;
        }while(p);

        stmt_table_tail->next = malloc(sizeof(SP_STMT_LST));
        p = stmt_table_tail->next;
        p->ifam = stmt_table_tail->ifam + 1;
    }

    p->name = strdup(elem);
    p->mag  = NULL;
    p->property = NULL;
    p->next = NULL;
    
    stmt_table_tail = p;
    return p;
}



/*! \brief Get the index of a beamline.
 *
 * \param bl Beamline name.
 * \return 0 if not defined if < 0, else return the index in beam_lines, 
 * (location of elements, key = elements - 1)
 * 
 * beam_lines has size of 2*n, for a n elements beam line. 2*i is name, 2*i+1 is definition (a string of all elements).
 */
int beamline_def(char* bl)
{
    size_t n = n_beam_lines;
    size_t i = 0;
    for (i = 0; i < n; i += 2) {
        if(str_case_cmp(bl, beam_lines[i]) == 0) return i + 1;
    }
    return 0;
}

/*! \brief Duplicate a beam line
 */
char* beamline_dup(char* bl)
{
    /* char* ret; */
    size_t i = beamline_def(bl);
    if ( i > 0 ) {
        return strdup(beam_lines[i]);
    }else {
        yyerror("Undefined beam line");
    }

    return NULL;
}

/*! \brief Add a beam line to global beam_lines[]
 */
void beamline_add(char* bl, char* elem)
{
    int i = 0;
    char **buf;
    /* fprintf(stderr, "N of Lines: %d\n", n_beam_lines); */

    /* check if exist */
    i = beamline_def(bl);
    if (i > 0) {
        free(beam_lines[i]);
        beam_lines[i] = strdup(elem);
    } else {
    /* add a new one */

        n_beam_lines += 2;
        buf = (char**) malloc(sizeof(char*)*n_beam_lines);
        for (i = 0; i < n_beam_lines - 2; ++i) {
            buf[i] = beam_lines[i];
        }
        free(beam_lines);

        beam_lines = buf;
        beam_lines[n_beam_lines-2] = strdup(bl);
        beam_lines[n_beam_lines-1] = strdup(elem);
    }
}

/*! \brief Reverse a string.
 */
static void reverse_str(char* s, int n)
{
    char c;
    char *p = s;
    int i = 0;
    for (i = 0; i < n/2; ++i) {
        c = *(p+i);
        *(p+i) = *(p+n-i-1);
        *(p+n-i-1) = c;
    }
}

/*! \brief Reverse a beam line.
 * (QF, DL, QD) -> (QD, DL, QF)
 */
char* beamline_reverse(char* bl)
{
    size_t i = 0, j =0;
    char *ret = strdup(bl);
    reverse_str(ret, strlen(ret));
    
    for (i = 0; i < strlen(ret); ++i) {
        if (*(ret+i) != ',') continue;
        reverse_str(ret+j, i-j);
        j = i+1;
    }

    reverse_str(ret+j, strlen(ret) - j);
    /* fprintf(stderr, "\n%s\n", bl);
       fprintf(stderr, "%s\n", ret); */
    return ret;
}


void free_beam_lines()
{
    int i = 0; 
    for (i = 0; i < n_beam_lines; ++i) {
        free(beam_lines[i]);
    }
    free(beam_lines);
    n_beam_lines = 0;
    beam_lines = NULL;
}


/* -------------------------------------------------- */

/* extern FILE* yyin; */

/*! \brief Print elements properties */
int print_elements(SP_STMT_LST* pelem)
{
    /* SP_STMT_LST *pelem = stmt_table; */
    SP_PRPT_LST *prpt;
    int i = 0;
    fprintf(stdout, "----------- elements -------------\n");
    while (pelem) { 
        fprintf(stdout, "%s: %s, ", pelem->name, pelem->mag);
        prpt = pelem->property;
        while (prpt) {
            if (prpt->str) {
                fprintf(stdout, " %s=\"%s\"", prpt->property, prpt->str);
            }else {
            fprintf(stdout, " %s= ( ", prpt->property);
            for (i=0; i < prpt->nval - 1; ++i)
                fprintf(stdout, "%f, ", prpt->val[i]);
            fprintf(stdout, "%f ) ", prpt->val[prpt->nval-1]);
            }
            prpt = prpt->next;
        }
        pelem = pelem->next;
        fprintf(stdout, "\n");
    }
}

int parse_lattice(const char* f)
{
    /* add mathematical function support */
    extern double sqrt(), exp(), log(), atan(), sin(), cos(), fabs();

    /* initialize the random generator */
    srand(0);

    if (!(yyin = fopen(f, "r"))) {
        fprintf(stderr, "ERROR: can  not open file %s\n", f);
        return -1;
    }


    sym_addfunc("sqrt", sqrt);
    sym_addfunc("exp", exp);
    sym_addfunc("log", log);
    sym_addfunc("log10", log10);

    sym_addfunc("sin", sin);
    sym_addfunc("cos", cos);
    sym_addfunc("tan", tan);
    sym_addfunc("asin", asin);
    sym_addfunc("acos", acos);
    sym_addfunc("atan", atan);
    sym_addfunc("cosh", cosh);
    sym_addfunc("sinh", sinh);
    sym_addfunc("tanh", tanh);

    sym_addfunc("abs", fabs);

    /* yypush_buffer_state(yycreate_buffer(yyin, 10240)); */

    while(!feof(yyin)) {
        yyparse();
        /* fprintf(stdout, "Returned from parser!\n"); */
    }
    fclose(yyin);

 #ifdef DEBUG
    fprintf(stdout, "\n#------------------------------------------------#\n");
    fprintf(stdout, "# DEBUG: printing symbol table\n");
    fprintf(stdout, "#------------------------------------------------#\n\n");
    print_symb_tab(symtab);
    fprintf(stdout, "#------------------------------------------------#\n");

    fprintf(stdout, "\n#------------------------------------------------#\n");
    fprintf(stdout, "# DEBUG: printing all statements\n");
    fprintf(stdout, "#------------------------------------------------#\n\n");
    print_elements(stmt_table);
    fprintf(stdout, "#------------------------------------------------#\n");
 #endif
}

int free_lattice()
{
    /* Free all the memory for lattice parsing */
    free_stmt_table();
    free_beam_lines();

    free_symb_tab(symtab); 
    symtab = NULL;

    free(lat_title); 
    lat_title = NULL;
}

/*! \brief Parse a lattice file
 *
 * \param f lattice file name.
 * \return 0/1
 */
int parse_lattice_flat(char* f, char* flat)
{
    int i = 0;
    int n_elem = 0;
    char* str;

    FILE* pf;
    parse_lattice(f);

    /* output every element */
    /* print_elements(stmt_table); */
    if(echo && get_statement_property_text("set", "output", &str)) {
        fprintf(stdout, "Output: %s %zd\n", str, strlen(str));
    }
    else fprintf(stdout, "Not found output\n");
    
    /*
    double d, *v;
    if (get_statement_property_double("set", "energy", &d)) 
        fprintf(stdout, "energy: %f\n", d);
    else fprintf(stdout, "Not found energy\n");

    if (get_sym_double_vec("C", &i, &v)) {
        fprintf(stdout, "C  %d\n", i);
        for (j = 0; j < i; ++j) fprintf(stdout, " %f", v[j]);
        fprintf(stdout, "\n");
        free(v);
    }
    */

    /* printf("res: %s\n", sp->str); */
    #ifdef USE_DB
    save_elements_db(stmt_table, str);
    #endif
    free(str);
    
    /* output every beam line */
    for (i = 0; i < n_beam_lines; i+=2) {
        char *p = beam_lines[i+1];
        while(*p) {
            if (*p == ',') ++n_elem;
            ++p;
        }
        ++n_elem;
        /* fprintf(stdout, "%s = %s; [%d]\n", beam_lines[i], beam_lines[i+1], n_elem); */
        n_elem = 0;
    }

    if ( n_beam_lines <= 0) {
        fprintf(stderr, "Error: No beamline specified in latticefile.\n");
        return 0;
    }

    pf = fopen(flat, "w");
    if (pf) {
        /* use the last beamlines as the "RING" */
        print_flat_lattice(pf, beam_lines[n_beam_lines-1]);
        fclose(pf);
    } else if (n_beam_lines > 0) {
        fprintf(stderr, "Error: can not open %s to write flat file.\n", flat);
        exit(-1);
    } 

    free_lattice();
    return 0;
}


/*! \brief get the path and root file name 
 *
 * allocated new mem, need free by caller.
 */
char *get_path_root(char *f)
{
    size_t n = 0; //
    char *r;
    char* s = f + strlen(f);
    while(s > f && *s != '.') {
        /* fprintf(stdout, "%c", *s); */
        --s;
    }
    /* fprintf(stdout, "\n"); */

    if (s == f) return NULL;
    n = s - f;

    r = malloc(sizeof(char)*(n+1));
    strncpy(r, f, n);
    r[n] = '\0';

    return r;
}


/*! \brief Parse a lattice file
 *
 * \param f lattice file name.
 * \return 0/1
 */
int parse_lattice_xml(char* lat, char* xml)
{
    int i = 0;
    int n_elem = 0;
    char *str;
    FILE *pf;
    int vmaj, vmin, vpatch, vrevision;
    SP_STMT_LST *pstmt;
    SP_PRPT_LST *prpt;
    
    /*
    str = get_path_root(f);
    xml = malloc(sizeof(char)*(strlen(str)+5));
    strcpy(xml, str);
    strcat(xml, ".xml");
    free(str);
    */

    if ((pf = fopen(xml, "wt")) == NULL) {
        fprintf(stderr, "[ERROR] Can not open XML file to write\n");
        return -1;
    };
    
    fprintf(stderr, "XML: %s\n", xml);

    get_parser_version(&vmaj, &vmin, &vpatch, &vrevision);

    if (parse_lattice(lat)) {
        // did not succeed
        fprintf(stderr, "Can not parse lattice\n");
        return -1;
    }

    fprintf(pf, "<?xml version=\"1.0\" ?>\n");
    fprintf(pf, "<latparser version=\"%d.%d.%d\" revision=\"%d\" build=\"%s\" />\n", 
            vmaj, vmin, vpatch, vrevision, __DATE__);
    /* output every element */
    /* print_elements(stmt_table); */
    fprintf(pf, "<lattice>\n");
    pstmt = get_statement_list();

    while(pstmt) {
        // loop over all statements
        fprintf(pf, "  <statement family=\"%s\" type=%d ifam=%d",
                pstmt->mag, pstmt->type, pstmt->ifam);
        prpt = pstmt->property;

        while (prpt) {
            fprintf(pf, " %s=\"", prpt->property);
            if (prpt->nval > 0) {
                for (i = 0; i < prpt->nval-1; ++i)
                    fprintf(pf, "%-.16E,", prpt->val[i]);
                fprintf(pf, "%-.16E\"", prpt->val[prpt->nval-1]);
            } else {
                fprintf(pf, "%s\"", prpt->str);
            }
            prpt = prpt->next;
        }
        fprintf(pf, ">%s</statement>\n", pstmt->name);
        pstmt = pstmt->next;
    }

    fprintf(pf, "</lattice>\n");

    free_lattice();
    fclose(pf);
    return 0;
}


/*! \brief Parse a lattice file
 *
 * \param f lattice file name.
 * \return 0/1
 */
int parse_lattice_ini(char* f)
{
    int i = 0;
    int n_elem = 0;
    char* str;
	char *ini;
    FILE* pf;
    int vmaj, vmin, vpatch, vrevision;
    SP_STMT_LST *pstmt;
    SP_PRPT_LST *prpt;
    
    str = get_path_root(f);

    ini = malloc(sizeof(char)*(strlen(str)+5));
    strcpy(ini, str);
    strcat(ini, ".ini");
    free(str);

    fprintf(stderr, "INI: %s\n", ini);

    get_parser_version(&vmaj, &vmin, &vpatch, &vrevision);

    if (parse_lattice(f)) {
        // did not succeed
        fprintf(stderr, "Can not parse lattice\n");
        return -1;
    }

    pf = fopen(ini, "w");
    if (pf == NULL) {
        fprintf(stderr, "Can not open %s\n", ini);
        return -1;
    }
    fprintf(pf, "[__info__]\n"
            "  version=%d.%d.%d\n  revision=%d\n  build=%s\n\n", 
            vmaj, vmin, vpatch, vrevision, __DATE__);
    /* output every element */
    /* print_elements(stmt_table); */
    pstmt = get_statement_list();

    while(pstmt) {
        // loop over all statements
        fprintf(pf, "[%s]\n", pstmt->name);
        fprintf(pf, "  magnet = %s\n", pstmt->mag);
        fprintf(pf, "  type = %d\n", pstmt->type);
        fprintf(pf, "  ifam = %d\n", pstmt->ifam);

        prpt = pstmt->property;

        while (prpt) {
            fprintf(pf, "  %s = ", prpt->property);
            if (prpt->nval > 0) {
                for (i = 0; i < prpt->nval-1; ++i)
                    fprintf(pf, "%f,", prpt->val[i]);
                fprintf(pf, "%f\n", prpt->val[prpt->nval-1]);
            } else {
                fprintf(pf, "%s\n", prpt->str);
            }
            /* fprintf(pf, "\n"); */
            prpt = prpt->next;
        }
        fprintf(pf, "\n");
        pstmt = pstmt->next;
    }

    fprintf(pf, "\n");
    fclose(pf);

    free_lattice();

    return 0;
}

/*! \brief Print the element into a flat machine file format.
 */
void print_flat_element(FILE* pf, SP_STMT_LST *elem)
{
    SP_STMT_LST *p;
    SP_PRPT_LST *prpt;
    int i, method = 4, N;
    double L = 0.0, lambda = 0.0;
    double xmin = -10., xmax = 10., ymin = -10., ymax = 10.;
    double rho = 0.0, curv = 0.0, t, t1 = 0.0, t2=0.0;
    double c0 = 2.99792458e8;

    double dx = 0.0, dy = 0.0, dtheta = 0.0;
    int n_mult_coeff = 0;

    double K = 0.0;
    double energy = 3.0;
    int h = 1300;
    /* type code */
    fprintf(pf, " %3d", elem->type);
    
    prpt = get_statement_property(elem, "method");
    method = prpt->nval ? prpt->val[0] : 0;
    
    prpt = get_statement_property(elem, "n");
    N = prpt->nval ? prpt->val[0] : 0;

    prpt = get_statement_property(elem, "l");
    L = prpt->nval ? prpt->val[0]:0.0;



    switch(elem->type) {
    case ELEMT_MARKER:
        /* integration method, num of integration steps. */
        fprintf(pf, " %3d %3d\n", method, N);
        /* apertures. */
        fprintf(pf, " %23.16e %23.16e %23.16e %23.16e\n",
	        xmin, xmax, ymin, ymax);
        break;
    case ELEMT_DRIFT:
        /* integration method, num of integration steps. */
        fprintf(pf, " %3d %3d\n", method, N);
        fprintf(pf, " %23.16e %23.16e %23.16e %23.16e\n",
                xmin, xmax, ymin, ymax);
        fprintf(pf, " %23.16e\n", L);
        break;
    case ELEMT_MULTIPOLE:
        /* integration method, num of integration steps. */
        fprintf(pf, " %3d %3d\n", method, N);
        /* apertures: xmin, xmax, ymin, ymax */
        fprintf(pf, " %23.16e %23.16e %23.16e %23.16e\n",
                xmin, xmax, ymin, ymax);
        /* h/v displacement, roll angle(design), roll angle(error). */
        fprintf(pf, " %23.16e %23.16e %23.16e %23.16e\n",
                dx, dy, dtheta, 0.0);

        if (str_case_cmp(elem->mag, "bend") == 0) {
            prpt = get_statement_property(elem, "t");
            curv = prpt->nval ? (prpt->val[0]*M_PI/180)/L : 0.0;
            prpt = get_statement_property(elem, "t1");
            t1 = prpt->nval ? prpt->val[0] : 0.0;
            prpt = get_statement_property(elem, "t2");
            t2 = prpt->nval ? prpt->val[0] : 0.0;
            fprintf(pf, " %23.16e %23.16e %23.16e %23.16e %23.16e\n",
                    L, curv, t1, t2, 0.0);
            /* num of nonzero multipole coeff. */
            fprintf(pf, "  %1d\n", 0);
            break;
        }
        prpt = get_statement_property(elem, "k");
        K = prpt->nval ? prpt->val[0] : 0.0;

        if (str_case_cmp(elem->mag, "quadrupole") == 0){
            i = 2;
            fprintf(pf, " %23.16e %23.16e %23.16e %23.16e %23.16e\n",
                    L, curv, 0.0, 0.0, 0.0);
            /* num of nonzero multipole coeff. */
            fprintf(pf, "  %1d %1d\n", 1, 2);
            fprintf(pf, "%3d %23.16e %23.16e\n", i, K, 0.0);
        }
        else if (str_case_cmp(elem->mag, "sextupole") == 0){
            i = 3;
            fprintf(pf, " %23.16e %23.16e %23.16e %23.16e %23.16e\n",
                    L, curv, 0.0, 0.0, 0.0);
            fprintf(pf, "  %1d %1d\n", 1, 3);
            fprintf(pf, "%3d %23.16e %23.16e\n", i, K, 0.0);
        }
        else {
            i = 0;
            fprintf(pf, "  %1d\n", 0);
        }
        break;
    case ELEMT_CAVITY: {
        double v, f;
        fprintf(pf, " %3d %3d\n", method, N);
        fprintf(pf, " %23.16e %23.16e %23.16e %23.16e\n",
                xmin, xmax, ymin, ymax);
        prpt = get_statement_property(elem, "frequency");
        f = prpt->nval ? prpt->val[0] : 0.0;
        prpt = get_statement_property(elem, "voltage");
        v = prpt->nval ? prpt->val[0] : 0.0;

        fprintf(pf, " %23.16e %23.16e %1d %23.16e\n",
                v/(1.0e9*energy), 2*M_PI*f/c0, h, 1e9*energy);
        break;}
    case ELEMT_THINKICK:
        if (method == 0) method = 4;
        fprintf(pf, " %3d %3d\n", method, N);
        fprintf(pf, " %23.16e %23.16e %23.16e %23.16e\n",
                xmin, xmax, ymin, ymax);
        fprintf(pf, " %23.16e %23.16e %23.16e\n",
                dx, dy, dtheta);
        fprintf(pf, "%3d\n", n_mult_coeff);
        break;
    case ELEMT_WIGGLER:
        fprintf(pf, " %3d %3d\n", method, N);
        fprintf(pf, " %23.16e %23.16e %23.16e %23.16e\n",
                xmin, xmax, ymin, ymax);
        prpt = get_statement_property(elem, "lambda");
        lambda = prpt->nval ? prpt->val[0] : 0.0;
        fprintf(pf, " %23.16e %23.16e\n", L, lambda);
        prpt = get_statement_property(elem, "harm");
        fprintf(pf, " %1zd\n", prpt->nval/6);

        for (i = 0; i < prpt->nval/6; ++i) {
            fprintf(pf, "%2d %23.16e %23.16e %23.16e %23.16e %23.16e\n",
                    (int)prpt->val[i*6], prpt->val[i*6+1],
                    prpt->val[i*6+2], prpt->val[i*6+3],
                    prpt->val[i*6+4], prpt->val[i*6+5]);
        }  
		  
        break;
    }
}


/*! \brief Print a beam line into flat machine file
 */
int print_flat_lattice(FILE* pf, char* bl)
{
    char *s = (char*) malloc(sizeof(char)*(strlen(bl)+1));
    /* fprintf(pf, "%s\n", bl); */
    int k, i, j, n;
    int n_elem = 0;
	char **elem;
	int *child;
	
    SP_STMT_LST *pelem;
    n = strlen(bl);

    elem = (char**) malloc(sizeof(char*)*n); /* element num in beam line*/
    child = (int*) malloc(sizeof(int)*n); /* kid number */


    for (k = i = 0; i < n + 1; ++i) {
        if (bl[i] != ',' && bl[i] != '\0') continue;
  
        strncpy(s, bl+k, i-k);
        s[i-k] = '\0';
        
        elem[n_elem] = strdup(s);
        child[n_elem%n] = 1;
        for (j = n_elem-1; j >= 0; --j) {
            if (str_case_cmp(elem[j], s) == 0) {
                child[n_elem%n] = child[j] + 1;
                break;
            }
        }
       
        
        pelem = stmt_table;
        while(pelem) {
            if (str_case_cmp(s, pelem->name) == 0) break;
            pelem = pelem->next;
        }

        fprintf(pf, "%-15s %4d %4d %4d\n", 
                s, pelem->ifam, child[n_elem], n_elem+1);

        print_flat_element(pf, pelem);

        ++n_elem;

        k = i + 1;
        ++i;
    }

    free(child);

    for (i = 0; i < n_elem; ++i) 
        free(elem[i]);
    free(elem);
    free(s);

    return 0;
}

SP_STMT_LST* get_statement_list() 
{
    return stmt_table;
}
