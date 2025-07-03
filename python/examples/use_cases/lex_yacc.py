# To install PLY:
#   pip install ply

import ply.lex as lex
import ply.yacc as yacc

#-------------------------------------------------------------------------------

# Lex.

def t_ID(t):
    r'[a-zA-Z_][a-zA-Z_0-9]*'
    t.type = reserved.get(t.value, 'ID')  # Check for reserved words
    return t

def t_NUMBER(t):
    r'\d+'
    t.value = int(t.value)
    return t

def t_newline(t):
    r'\n+'
    t.lexer.lineno += len(t.value)

def t_error(t):
    print(f"Illegal character '{t.value[0]}'")
    t.lexer.skip(1)

#-------------------------------------------------------------------------------

# Yacc.

def p_statement_expr(p):
    'statement : expression'
    print(p[1])

def p_statement_assign(p):
    'statement : ID EQUALS expression'
    variables[p[1]] = p[3]

def p_statement_print(p):
    'statement : PRINT ID'
    if p[2] in variables:
        print(variables[p[2]])
    else:
        print(f"Undefined variable '{p[2]}'")

def p_expression_binop(p):
    '''expression : expression PLUS term
                  | expression MINUS term'''
    if p[2] == '+': p[0] = p[1] + p[3]
    elif p[2] == '-': p[0] = p[1] - p[3]

def p_expression_term(p):
    'expression : term'
    p[0] = p[1]

def p_term_binop(p):
    '''term : term TIMES factor
            | term DIVIDE factor'''
    if p[2] == '*': p[0] = p[1] * p[3]
    elif p[2] == '/': p[0] = p[1] / p[3]

def p_term_factor(p):
    'term : factor'
    p[0] = p[1]

def p_factor_number(p):
    'factor : NUMBER'
    p[0] = p[1]

def p_factor_var(p):
    'factor : ID'
    if p[1] in variables:
        p[0] = variables[p[1]]
    else:
        print(f"Undefined variable '{p[1]}'")
        p[0] = 0

def p_factor_group(p):
    'factor : LPAREN expression RPAREN'
    p[0] = p[2]

def p_error(p):
    print("Syntax error!")

#-------------------------------------------------------------------------------

# Lex.
# List of token names
tokens = (
    'NUMBER',
    'ID',
    'PLUS', 'MINUS', 'TIMES', 'DIVIDE',
    'EQUALS',
    'LPAREN', 'RPAREN',
    'PRINT',
)

# Regular expression rules
t_PLUS    = r'\+'
t_MINUS   = r'-'
t_TIMES   = r'\*'
t_DIVIDE  = r'/'
t_EQUALS  = r'='
t_LPAREN  = r'\('
t_RPAREN  = r'\)'
t_ignore  = ' \t'

# Reserved words
reserved = {
    'print': 'PRINT',
}


# Yacc.
# Symbol table (for storing variables)
variables = {}

lexer = lex.lex()
parser = yacc.yacc()

while True:
    try:
        s = input('>>> ')
    except EOFError:
        break
    parser.parse(s)

# Example:
# x = 4 + 5
# y = x * 2
# print y

