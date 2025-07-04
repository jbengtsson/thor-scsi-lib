
import sys
from pathlib import Path

from sly import Lexer, Parser


# Dummy placeholders for external functions and types:
# Replace these with your actual implementations

def glps_string_alloc(text):
    return text

def glps_error(scanner, context, msg, *args):
    if args:
        msg = msg % args
    raise SyntaxError(msg)

def glps_assign(context, key, value):
    print(f"[assign] {key} = {value}")
    context['assignments'][key] = value

def glps_add_element(context, cat, name, props):
    print(f"[element] {cat}:{name} with props {props}")
    context['elements'].append(
        {'category': cat, 'name': name, 'properties': props})

def glps_call1(context, func, arg):
    print(f"[func] {func}({arg})")
    context['functions'].append({'func': func, 'arg': arg})

def glps_command(context, cmd):
    print(f"[command] {cmd}")
    context['commands'].append(cmd)

def glps_expr_number(val): return {'type': 'number', 'value': val}
def glps_expr_string(val): return {'type': 'string', 'value': val}
def glps_expr_var(val): return {'type': 'var', 'name': val}
def glps_expr_vector(vec): return {'type': 'vector', 'items': vec}

def glps_add_value(context, kind_func, value):
    return kind_func(value)

def glps_add_op(context, op, arity, args):
    return {'type': 'op', 'op': op, 'args': args}


class GLPSLexer(Lexer):
    tokens = { 'IDENT', 'NUM', 'STR' }
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
        raw = t.value[1:-1]  # strip the quotes
        try:
            # Decode escape sequences manually
            t.value = glps_string_alloc(
                bytes(raw, "utf-8").decode("unicode_escape"))
        except Exception as e:
            glps_error(None, None, f"Invalid string escape: {e}")
        return t

    @_(r'\n+')
    def ignore_newline(self, t):
        self.lineno += t.value.count('\n')
        self.line_start = t.index + 1

    def error(self, t):
        line = self.lineno
        col = t.index - self.line_start if hasattr(t, 'index') else 0
        msg = f"Invalid character {t.value[0]!r} at line {line}, column {col}"
        # Skip the invalid character by advancing the lexer position manually
        # Skip invalid character to avoid infinte loops.
        self.index = t.index + 1
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

    @_('assignment', 'element', 'func', 'command')
    def entry(self, p):
        return p[0]

    @_('IDENT "=" expr ";"')
    def assignment(self, p):
        glps_assign(self.ctxt, p.IDENT, p.expr)
        return ('assign', p.IDENT, p.expr)

    # No optional properties.
    @_('IDENT ":" IDENT ";"')
    def element(self, p):
        glps_add_element(self.ctxt, p.IDENT0, p.IDENT1, [])
        return ('element', p.IDENT0, p.IDENT1, [])

    # Leading comma required for optional properties.
    @_('IDENT ":" IDENT "," property_list ";"')
    def element(self, p):
        glps_add_element(self.ctxt, p.IDENT0, p.IDENT1, p.property_list)
        return ('element', p.IDENT0, p.IDENT1, p.property_list)

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

    @_('IDENT "(" expr ")" ";"')
    def func(self, p):
        glps_call1(self.ctxt, p.IDENT, p.expr)
        return ('func', p.IDENT, p.expr)

    @_('IDENT ";"')
    def command(self, p):
        glps_command(self.ctxt, p.IDENT)
        return ('command', p.IDENT)

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
            line = getattr(p, 'lineno', None)
            col = None
            if hasattr(p, 'lexpos') and self.lexer \
               and hasattr(self.lexer, 'line_start'):
                col = p.lexpos - self.lexer.line_start + 1  # +1 for 1-based column
            msg = f"Syntax error at token '{value}' (type {p.type})"
            if line is not None:
                msg += f" on line {line}"
            if col is not None:
                msg += f", column {col}"
        else:
            msg = "Syntax error: unexpected end of input"
        print(msg)
        raise SyntaxError(msg)


if __name__ == "__main__":

    context = {
        'assignments': {},
        'elements': [],
        'functions': [],
        'commands': []
    }

    # Assign the quote lexer subclass to the main lexer class attribute

    lexer = GLPSLexer()
    parser = GLPSParser(ctxt=context, lexer=lexer)

    home_dir = Path.home() / "Nextcloud" / "thor_scsi" / "JB" / "MAX_IV"

    if len(sys.argv) < 2:
        print("Usage: script.py <filename>")
        sys.exit(1)

    file_name = Path(home_dir) / (sys.argv[1]+".lat")

    try:
        with open(file_name, "r") as file:
            text = file.read()
    except FileNotFoundError:
        print(f"File {file_name} not found")
        sys.exit(1)

    try:
        result = parser.parse(lexer.tokenize(text))
        if False:
            print("\nParsed context\n\nassignments:")
            print(context['assignments'])
            print("\nelements:")
            print(context['elements'])
            print("\nfunctions:")
            print(context['functions'])
            print("\ncommands:")
            print(context['commands'])
    except SyntaxError as e:
        print(f"Syntax error: {e}")
