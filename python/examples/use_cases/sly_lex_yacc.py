# To install SLY:
#   pip install sly
from sly import Lexer, Parser


class CalcLexer(Lexer):
    tokens = {
        NAME, NUMBER, PLUS, MINUS, TIMES, DIVIDE, ASSIGN, LPAREN, RPAREN
    }

    ignore = ' \t'

    # Token regex.
    PLUS        = r'\+'
    MINUS       = r'-'
    TIMES       = r'\*'
    DIVIDE      = r'/'
    ASSIGN      = r'='
    LPAREN      = r'\('
    RPAREN      = r'\)'

    @_(r'[a-zA-Z_][a-zA-Z0-9_]*')
    def NAME(self, t):
        return t

    @_(r'\d+')
    def NUMBER(self, t):
        t.value = int(t.value)
        return t

    def error(self, t):
        print(f"Illegal character '{t.value[0]}'")
        self.index += 1

class CalcParser(Parser):
    tokens = CalcLexer.tokens

    # Operator precedence and associativity
    precedence = (
        ('left', PLUS, MINUS),
        ('left', TIMES, DIVIDE),
    )

    def __init__(self):
        self.env = {}

    @_('NAME ASSIGN expr')
    def statement(self, p):
        self.env[p.NAME] = p.expr
        return p.expr

    @_('expr')
    def statement(self, p):
        return p.expr

    @_('expr PLUS expr')
    @_('expr MINUS expr')
    @_('expr TIMES expr')
    @_('expr DIVIDE expr')
    def expr(self, p):
        if p[1] == '+':
            return p.expr0 + p.expr1
        elif p[1] == '-':
            return p.expr0 - p.expr1
        elif p[1] == '*':
            return p.expr0 * p.expr1
        elif p[1] == '/':
            return p.expr0 / p.expr1

    @_('LPAREN expr RPAREN')
    def expr(self, p):
        return p.expr

    @_('NUMBER')
    def expr(self, p):
        return p.NUMBER

    @_('NAME')
    def expr(self, p):
        try:
            return self.env[p.NAME]
        except KeyError:
            print(f"Undefined variable '{p.NAME}'")
            return 0


if False:
    # Test lexical analyser.
    lexer = CalcLexer()
    for tok in lexer.tokenize("a=3+4*5"):
        print(tok)

# Calculator.
if __name__ == '__main__':
    lexer = CalcLexer()
    parser = CalcParser()

    while True:
        try:
            text = input('calc > ')
            if text.strip() == '':
                continue
        except EOFError:
            break
        result = parser.parse(lexer.tokenize(text))
        if result is not None:
            print(result)

# Example:
#   calc > x = 4
#   4
#   calc > y = x * 2
#   8
#   calc > y + 6
#   14
#   calc > z
