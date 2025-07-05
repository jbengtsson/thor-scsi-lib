
import sys
from pathlib import Path

from sly import Lexer, Parser


# Dummy placeholders for external functions and types:
# Replace these with your actual implementations

debug = False

def glps_string_alloc(text):
    return text

def glps_error(scanner, ctxt, msg, *args):
    if args:
        msg = msg % args
    raise SyntaxError(msg)

def glps_assign(ctxt, key, value):
    if debug:
        print(f"[assign] {key} = {value}")
    ctxt['assignments'][key] = value

def glps_add_element(ctxt, name, type, props):
    if props is None:
        props = []
    if not isinstance(props, list):
        raise TypeError(f"props should be a list, got {type(props).__name__}")
    if debug:
        print(f"[element] {name}:{type} with props {props}")
    ctxt['elements'].append(
        {'name': name, 'type': type, 'properties': props})

def glps_add_line(ctxt, name, key2, line_list):
    if debug:
        print(f"[line] {name}:{key2} = {line_list}")
    ctxt['lines'][name] = line_list

def glps_call1(ctxt, func, arg):
    if debug:
        print(f"[func] {func}({arg})")
    ctxt['functions'].append({'func': func, 'arg': arg})

def glps_command(ctxt, cmd, arg):
    if debug:
        print(f"[command] {cmd}: {arg}")
    ctxt['commands'][cmd] = arg

def glps_append_expr(ctxt, expr_list, expr):
    if expr_list is None:
        expr_list = []
    expr_list.insert(0, expr)  # emulate right-recursion order
    return expr_list

def glps_add_value(ctxt, type_, val):
    return (type_, val)

def glps_add_op(ctxt, op, arity, args):
    return {'type': 'op', 'op': op, 'args': args}

# Types for demonstration:
glps_expr_number = 'number'
glps_expr_vector = 'vector'
glps_expr_string = 'string'
glps_expr_var    = 'elem'


class GLPSLexer(Lexer):
    tokens = { 'IDENT', 'NUM', 'STR', 'USE' }
    literals = { '=', ':', ';', '(', ')', '[', ']', ',', '+', '-', '*', '/' }
    # Ignore white space.
    ignore = ' \t'
    # Ignore everything from "#" to end-of-line.
    ignore_comment = r'\#[^\n]*'

    def __init__(self):
        super().__init__()
        self.line_start = 0  # Track start index of current line

    @_(r'[A-Za-z_][A-Za-z0-9_]*')
    def IDENT(self, t):
        if t.value.upper() == "USE":
            t.type = 'USE'
        t.value = glps_string_alloc(t.value)
        return t

    @_(r'[0-9]+(\.[0-9]*)?([eE][+-]?[0-9]+)?')
    def NUM(self, t):
        try:
            t.value = float(t.value)
        except ValueError:
            glps_error(None, None, "Invalid number: %s", t.value)
        return t

    @_(r'\"([^\"\\\n\r]|\\.)*\"')
    def STR(self, t):
        raw = t.value[1:-1]  # strip the quotes.
        try:
            # Decode escape sequences manually.
            t.value = glps_string_alloc(
                bytes(raw, "utf-8").decode("unicode_escape"))
        except Exception as e:
            glps_error(None, None, f"Invalid string escape: {e}")
        return t

    @_(r'\n+')
    def ignore_newline(self, t):
        self.lineno += t.value.count('\n')
        # self.line_start = self.index
        self.line_start = t.index + len(t.value)

    def error(self, t):
        line = self.lineno
        col = t.lexpos - self.line_start
        msg = f"Invalid character {t.value[0]!r} at line {line}, column {col}"
        # Skip invalid character to avoid infinte loops.
        self.index = t.lexpos + 1
        glps_error(None, None, msg)

        
class GLPSParser(Parser):
    tokens = GLPSLexer.tokens
    precedence = (
        ('left', '+', '-'),
        ('left', '*', '/'),
        ('right', 'NEG')
    )

    def __init__(self, ctxt, lexer):
        self.ctxt = ctxt
        self.lexer = lexer  # Save lexer reference for error reporting.

    @_('entries')
    def file(self, p):
        return {'context': self.ctxt, 'entries': p.entries}

    @_('entry entries')
    def entries(self, p):
        return [p.entry] + p.entries

    @_('entry')
    def entries(self, p):
        return [p.entry]

    @_('assignment', 'element', 'line', 'func', 'command')
    def entry(self, p):
        return p[0]

    @_('IDENT "=" expr ";"')
    def assignment(self, p):
        glps_assign(self.ctxt, p.IDENT, p.expr)
        return ('assign', p.IDENT, p.expr)

    # Leading comma required for optional properties.
    @_('IDENT ":" IDENT property_opt ";"')
    def element(self, p):
        props = p.property_opt if p.property_opt is not None else []
        glps_add_element(self.ctxt, p.IDENT0, p.IDENT1, props)
        return ('element', p.IDENT0, p.IDENT1, props)

    @_('"," property_list')
    def property_opt(self, p):
        return p.property_list

    @_('')
    def property_opt(self, p):
        return []

    # Property_list is a comma separated list with one or more properties.
    @_('property')
    def property_list(self, p):
        return [p.property]

    @_('property "," property_list')
    def property_list(self, p):
        return [p.property] + p.property_list

    @_('IDENT "=" expr')
    def property(self, p):
        return {'key': p.IDENT, 'value': p.expr}

    @_('IDENT ":" IDENT "=" "(" line_list ")" ";"')
    def line(self, p):
        glps_add_line(self.ctxt, p.IDENT0, p.IDENT1, p.line_list)
        return (p.IDENT0, p.IDENT1, p.line_list)

    @_('IDENT "(" expr ")" ";"')
    def func(self, p):
        glps_call1(self.ctxt, p.IDENT, p.expr)
        return ('func', p.IDENT, p.expr)

    @_('USE ":" IDENT ";"')
    def command(self, p):
        glps_command(self.ctxt, p.USE, p.IDENT)
        return ('command', p.USE, p.IDENT)

    # No IDENT without an argument.
    # @_('IDENT ";"')
    # def command(self, p):
    #     glps_command(self.ctxt, p.IDENT, '')
    #     return ('command', p.IDENT)

    @_('')
    def line_list(self, p):
        return None

    @_('expr')
    def line_list(self, p):
        return glps_append_expr(self.ctxt, None, p.expr)

    @_('expr "," line_list')
    def line_list(self, p):
        return glps_append_expr(self.ctxt, p.line_list, p.expr)

    @_('NUM')
    def expr(self, p):
        return glps_add_value(self.ctxt, glps_expr_number, p.NUM)

    @_('STR')
    def expr(self, p):
        return glps_add_value(self.ctxt, glps_expr_string, p.STR)

    @_('IDENT')
    def expr(self, p):
        return glps_add_value(self.ctxt, glps_expr_var, p.IDENT)

    @_('expr "+" expr',
       'expr "-" expr',
       'expr "*" expr',
       'expr "/" expr')
    def expr(self, p):
        return glps_add_op(self.ctxt, p[1], 2, [p[0], p[2]])

    @_('"-" expr %prec NEG')
    def expr(self, p):
        return glps_add_op(self.ctxt, '-', 1, [p.expr])

    @_('"(" expr ")"')
    def expr(self, p):
        return p.expr

    @_('"[" expr_list "]"')
    def expr(self, p):
        return glps_add_value(self.ctxt, glps_expr_vector, p.expr_list)

    @_('expr')
    def expr_list(self, p):
        return [p.expr]

    @_('expr "," expr_list')
    def expr_list(self, p):
        return [p.expr] + p.expr_list

    def error(self, p):
        if p:
            value = getattr(p, 'value', None)
            line = getattr(p, 'lineno', 'unknown')
            col = '?'

            # Try to get column number if possible
            if hasattr(p, 'lexpos') and hasattr(self.lexer, 'line_start'):
                try:
                    col = p.lexpos - self.lexer.line_start + 1
                except Exception:
                    col = '?'

            # Heuristic: check if it's a missing semicolon before IDENT
            if p.type == 'IDENT' and isinstance(value, str):
                msg = f"Syntax error: possible missing semicolon before " \
                    f"'{value}'"
            else:
                msg = f"Syntax error at token '{value}' (type {p.type})"

            msg += f" on line {line}, column {col}"
        else:
            msg = "Syntax error: unexpected end of input"

        print(msg)
        raise SyntaxError(msg)



def prt_assigns(assign):
    print("\nassignments:")
    for a in assign:
        # print(f"  {a:15s} {context['assignments'][a][1]:10.3e}")
        print(f"  {a:15s} {context['assignments'][a]}")

def prt_elem(elem):
    print(f"  {elem['name']:15s} {elem['type']:10s}", end="")
    if len(elem['properties']) != 0:
        print(f" {elem['properties'][0]['key']:10s}", end="")
        if isinstance(elem['properties'][0]['value'], tuple):
            print(f" {elem['properties'][0]['value'][1]:10.3e}")
        else:
            print()
    else:
        print()

def prt_elems(elem):
    print("\nelements:")
    for e in elem:
        prt_elem(e)

def prt_lines(line):
    n_prt = 5
    print("\nlines:")
    for l in line:
        print(f"{l:15s}\n   ", end="")
        for k, e in enumerate(line[l]):
            if k < n_prt:
                # Check for reverse elements.
                if isinstance(e, tuple):
                    print(f" {e[1]}", end="")
                else:
                    print(f" -()", end="")
            elif k == n_prt:
                print("\n    ...\n   ", end="")
            elif k >= len(line[l])-n_prt:
                # Check for reverse elements.
                if isinstance(e, tuple):
                    print(f" {e[1]}", end="")
                else:
                    print(f" -()", end="")
        print()

def prt_ctxt(context):
    prt_assigns(context['assignments'])
    prt_elems(context['elements'])
    prt_lines(context['lines'])
    print("\nfunctions:")
    print(context['functions'])
    print("\ncommands:")
    print(context['commands'])



if __name__ == "__main__":

    context = {
        'assignments': {},
        'elements':    [],
        'lines':       {},
        'functions':   [],
        'commands':    {}
    }

    # Assign the quote lexer subclass to the main lexer class attribute

    lexer = GLPSLexer()
    parser = GLPSParser(ctxt=context, lexer=lexer)

    home_dir = Path.home() / "Nextcloud" / "thor_scsi" / "JB" / "MAX_IV"

    if len(sys.argv) != 2:
        print("Usage: script.py <file name without extension>")
        sys.exit(1)

    file_name = Path(home_dir) / (sys.argv[1]+".lat")
    if not file_name.exists():
        print(f"File not found: {file_name}")
        sys.exit(1)

    try:
        with open(file_name, "r") as file:
            text = file.read()
    except FileNotFoundError:
        print(f"File {file_name} not found")
        sys.exit(1)

    try:
        result = parser.parse(lexer.tokenize(text))
        if False:
            print("\nAssignments:")
            for a in context['assignments']:
                print(f"{a:15s} {context['assignments'][a]}")
                print("\nLines:")
                for l in context['lines']:
                    print(f"\n {l:15s} {context['lines'][l]}")
        if not False:
            prt_ctxt(context)
    except SyntaxError as e:
        print(f"Syntax error: {e}")
