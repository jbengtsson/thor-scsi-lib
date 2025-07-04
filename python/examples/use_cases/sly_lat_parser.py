
import os
import sys
from sly import Lexer, Parser


# Dummy placeholders for external functions and types:
# Replace these with your actual implementations

def glps_string_alloc(text):
    return text

def glps_error(scanner, context, msg, *args):
    raise SyntaxError(msg % args)

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

def glps_append_kv(context, kvlist, kv):
    kvlist.insert(0, kv)
    return kvlist


class GLPSLexer(Lexer):
    tokens = { IDENT, NUM, STR }
    literals = { '=', ':', ';', '(', ')', '[', ']', ',', '+', '-', '*', '/' }
    ignore = ' \t'
    ignore_comment = r'\#[^\n]*'

    states = { ('quote', 'exclusive') }

    quote = None  # assigned after QuoteLexer

    def __init__(self):
        self._quoted_str = ""
        self.quote = QuoteLexer()

    @_(r'[A-Za-z_][A-Za-z0-9_:]*')
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

    @_(r'\n')
    def ignore_newline(self, t):
        self.lineno += 1

    @_(r'\"')
    def quote_start(self, t):
        self._quoted_str = ""
        self.begin('quote')

    def error(self, t):
        pos = t.index
        line_start = self.text.rfind('\n', 0, pos) + 1
        col = pos - line_start + 1
        msg = f"Invalid character {t.value[0]!r} at line {t.lineno}, column {col}"
        glps_error(None, None, msg)


class QuoteLexer(Lexer):
    tokens = { 'STR' }
    ignore = ''

    @_(r'\"')
    def quote_STR(self, t):
        self.begin('INITIAL')
        t.type = 'STR'
        t.value = glps_string_alloc(self._quoted_str)
        return t

    @_(r'[^\"\n\r]+')
    def quote_text(self, t):
        self._quoted_str += t.value

    @_(r'[\n\r]')
    def quote_error(self, t):
        glps_error(None, None, f"Unterminated string at line {t.lineno}")

    @_(r'.')
    def quote_error_general(self, t):
        glps_error(None, None, f"Invalid character in string: {t.value[0]!r}")


class GLPSParser(Parser):
    tokens = GLPSLexer.tokens
    precedence = (
        ('left', '+', '-'),
        ('left', '*', '/'),
        ('right', 'NEG')
    )

    def __init__(self, ctxt):
        self.ctxt = ctxt

    @_('')
    def file(self, p):
        pass

    @_('entry file')
    def file(self, p):
        return None

    @_('assignment', 'element', 'func', 'command')
    def entry(self, p): pass

    @_('IDENT "=" expr ";"')
    def assignment(self, p):
        glps_assign(self.ctxt, p.IDENT, p.expr)

    @_('IDENT ":" IDENT properties ";"')
    def element(self, p):
        glps_add_element(self.ctxt, p.IDENT0, p.IDENT1, p.properties)

    @_('IDENT "(" expr ")" ";"')
    def func(self, p):
        glps_call1(self.ctxt, p.IDENT, p.expr)

    @_('IDENT ";"')
    def command(self, p):
        glps_command(self.ctxt, p.IDENT)

    @_('')
    def properties(self, p):
        return []

    @_('"," property properties')
    def properties(self, p):
        return glps_append_kv(self.ctxt, p.properties, p.property)

    @_('IDENT "=" expr')
    def property(self, p):
        return {'key': p.IDENT, 'value': p.expr}

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

    @_('')
    def expr_list(self, p): return []

    @_('expr')
    def expr_list(self, p): return [p.expr]

    @_('expr "," expr_list')
    def expr_list(self, p): return [p.expr] + p.expr_list

    def error(self, p):
        if p:
            value = getattr(p, 'value', None)
            line = getattr(p, 'lineno', None)
            msg = f"Syntax error at token '{value}' (type {p.type})"
            if line is not None:
                msg += f" on line {line}"
            print(msg)
        else:
            print("Syntax error: unexpected end of input")
        raise SyntaxError("Parsing failed")


if __name__ == "__main__":

    context = {
        'assignments': {},
        'elements': [],
        'functions': [],
        'commands': []
    }

    # Assign the quote lexer subclass to the main lexer class attribute

    GLPSLexer.quote = QuoteLexer

    lexer = GLPSLexer()
    parser = GLPSParser(ctxt=context)

    home_dir = os.path.join(
        os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")

    if len(sys.argv) < 2:
        print("Usage: script.py <filename>")
        sys.exit(1)

    file_name = os.path.join(home_dir, sys.argv[1]+".lat")

    try:
        with open(file_name, "r") as file:
            text = file.read()
    except FileNotFoundError:
        print(f"File {file_name} not found")
        sys.exit(1)

    try:
        result = parser.parse(lexer.tokenize(text))
        print("\nParsed context\n\nassignments:")
        print(context['assignments'])
        print("\nelements:")
        print(context['elements'])
        print("\nfunctions:")
        print(context['functions'])
        print("\ncommands:")
        print(context['commands'])
        print("\nParse Tree:")
        print(result)
    except SyntaxError as e:
        print(f"Syntax error: {e}")
