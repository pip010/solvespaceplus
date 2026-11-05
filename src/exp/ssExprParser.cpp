#include "ssExpParser.h"

#include "ssExp.h"

namespace SolveSpace
{

    class ExprParser;

Expr *Expr::Parse(const std::string &input, std::string *error) {
return ExprParser::Parse(input, error);
}

Expr *Expr::From(const std::string &input, bool popUpError) {
    std::string error;
    Expr *e = ExprParser::Parse(input, &error);
    if(!e) {
        dbp("Parse/lex error: %s", error.c_str());
        if(popUpError) {
            Error("Not a valid number or expression: '%s'.\n%s.",
                  input.c_str(), error.c_str());
        }
    }
    return e;
}


ExprParser::Token ExprParser::Token::From(TokenType type, Expr *expr) {
    Token t;
    t.type = type;
    t.expr = expr;
    return t;
}

ExprParser::Token ExprParser::Token::From(TokenType type, Expr::Op op) {
    Token t;
    t.type = type;
    t.expr = AllocExpr();
    t.expr->op = op;
    return t;
}

char ExprParser::ReadChar() {
    return *it++;
}

char ExprParser::PeekChar() {
    if(it == end) {
        return '\0';
    } else {
        return *it;
    }
}

std::string ExprParser::ReadWord() {
    std::string s;

    while(char c = PeekChar()) {
        if(!isalnum(c)) break;
        s.push_back(ReadChar());
    }

    return s;
}

void ExprParser::SkipSpace() {
    while(char c = PeekChar()) {
        if(!isspace(c)) break;
        ReadChar();
    }
}

ExprParser::Token ExprParser::LexNumber(std::string *error) {
    std::string s;

    while(char c = PeekChar()) {
        if(!((c >= '0' && c <= '9') || c == 'e' || c == 'E' || c == '.' || c == '_')) break;
        if(c == '_') {
            ReadChar();
            continue;
        }
        s.push_back(ReadChar());
    }

    char *endptr;
    double d = strtod(s.c_str(), &endptr);

    Token t = Token::From();
    if(endptr == s.c_str() + s.size()) {
        t = Token::From(TokenType::OPERAND, Expr::Op::CONSTANT);
        t.expr->v = d;
    } else {
        *error = "'" + s + "' is not a valid number";
    }
    return t;
}

ExprParser::Token ExprParser::Lex(std::string *error) {
    SkipSpace();

    Token t = Token::From();
    char c = PeekChar();
    if(isupper(c)) {
        std::string n = ReadWord();
        t = Token::From(TokenType::OPERAND, Expr::Op::VARIABLE);
    } else if(isalpha(c)) {
        std::string s = ReadWord();
        if(s == "sqrt") {
            t = Token::From(TokenType::UNARY_OP, Expr::Op::SQRT);
        } else if(s == "square") {
            t = Token::From(TokenType::UNARY_OP, Expr::Op::SQUARE);
        } else if(s == "sin") {
            t = Token::From(TokenType::UNARY_OP, Expr::Op::SIN);
        } else if(s == "cos") {
            t = Token::From(TokenType::UNARY_OP, Expr::Op::COS);
        } else if(s == "asin") {
            t = Token::From(TokenType::UNARY_OP, Expr::Op::ASIN);
        } else if(s == "acos") {
            t = Token::From(TokenType::UNARY_OP, Expr::Op::ACOS);
        } else if(s == "pi") {
            t = Token::From(TokenType::OPERAND, Expr::Op::CONSTANT);
            t.expr->v = PI;
        } else {
            *error = "'" + s + "' is not a valid variable, function or constant";
        }
    } else if(isdigit(c) || c == '.') {
        return LexNumber(error);
    } else if(ispunct(c)) {
        ReadChar();
        if(c == '+') {
            t = Token::From(TokenType::BINARY_OP, Expr::Op::PLUS);
        } else if(c == '-') {
            t = Token::From(TokenType::BINARY_OP, Expr::Op::MINUS);
        } else if(c == '*') {
            t = Token::From(TokenType::BINARY_OP, Expr::Op::TIMES);
        } else if(c == '/') {
            t = Token::From(TokenType::BINARY_OP, Expr::Op::DIV);
        } else if(c == '(') {
            t = Token::From(TokenType::PAREN_LEFT);
        } else if(c == ')') {
            t = Token::From(TokenType::PAREN_RIGHT);
        } else {
            *error = "'" + std::string(1, c) + "' is not a valid operator";
        }
    } else if(c == '\0') {
        t = Token::From(TokenType::END);
    } else {
        *error = "Unexpected character '" + std::string(1, c) + "'";
    }

    return t;
}

ExprParser::Token ExprParser::PopOperand(std::string *error) {
    Token t = Token::From();
    if(stack.empty() || stack.back().type != TokenType::OPERAND) {
        *error = "Expected an operand";
    } else {
        t = stack.back();
        stack.pop_back();
    }
    return t;
}

ExprParser::Token ExprParser::PopOperator(std::string *error) {
    Token t = Token::From();
    if(stack.empty() || (stack.back().type != TokenType::UNARY_OP &&
                         stack.back().type != TokenType::BINARY_OP)) {
        *error = "Expected an operator";
    } else {
        t = stack.back();
        stack.pop_back();
    }
    return t;
}

int ExprParser::Precedence(Token t) {
    ssassert(t.type == TokenType::BINARY_OP ||
             t.type == TokenType::UNARY_OP ||
             t.type == TokenType::OPERAND,
             "Unexpected token type");

    if(t.type == TokenType::UNARY_OP) {
        return 30;
    } else if(t.expr->op == Expr::Op::TIMES ||
              t.expr->op == Expr::Op::DIV) {
        return 20;
    } else if(t.expr->op == Expr::Op::PLUS ||
              t.expr->op == Expr::Op::MINUS) {
        return 10;
    } else if(t.type == TokenType::OPERAND) {
        return 0;
    } else ssassert(false, "Unexpected operator");
}

bool ExprParser::Reduce(std::string *error) {
    Token a = PopOperand(error);
    if(a.IsError()) return false;

    Token op = PopOperator(error);
    if(op.IsError()) return false;

    Token r = Token::From(TokenType::OPERAND);
    switch(op.type) {
        case TokenType::BINARY_OP: {
            Token b = PopOperand(error);
            if(b.IsError()) return false;
            r.expr = b.expr->AnyOp(op.expr->op, a.expr);
            break;
        }

        case TokenType::UNARY_OP: {
            Expr *e = a.expr;
            switch(op.expr->op) {
                case Expr::Op::NEGATE: e = e->Negate(); break;
                case Expr::Op::SQRT:   e = e->Sqrt(); break;
                case Expr::Op::SQUARE: e = e->Times(e); break;
                case Expr::Op::SIN:    e = e->Times(Expr::From(PI/180))->Sin(); break;
                case Expr::Op::COS:    e = e->Times(Expr::From(PI/180))->Cos(); break;
                case Expr::Op::ASIN:   e = e->ASin()->Times(Expr::From(180/PI)); break;
                case Expr::Op::ACOS:   e = e->ACos()->Times(Expr::From(180/PI)); break;
                default: ssassert(false, "Unexpected unary operator");
            }
            r.expr = e;
            break;
        }

        default: ssassert(false, "Unexpected operator");
    }
    stack.push_back(r);

    return true;
}

bool ExprParser::Parse(std::string *error, size_t reduceUntil) {
    while(true) {
        Token t = Lex(error);
        switch(t.type) {
            case TokenType::ERROR:
                return false;

            case TokenType::END:
            case TokenType::PAREN_RIGHT:
                while(stack.size() > 1 + reduceUntil) {
                    if(!Reduce(error)) return false;
                }

                if(t.type == TokenType::PAREN_RIGHT) {
                    stack.push_back(t);
                }
                return true;

            case TokenType::PAREN_LEFT: {
                // sub-expression
                if(!Parse(error, /*reduceUntil=*/stack.size())) return false;

                if(stack.empty() || stack.back().type != TokenType::PAREN_RIGHT) {
                    *error = "Expected ')'";
                    return false;
                }
                stack.pop_back();
                break;
            }

            case TokenType::BINARY_OP:
                if((stack.size() > reduceUntil && stack.back().type != TokenType::OPERAND) ||
                   stack.size() == reduceUntil) {
                    if(t.expr->op == Expr::Op::MINUS) {
                        t.type = TokenType::UNARY_OP;
                        t.expr->op = Expr::Op::NEGATE;
                        stack.push_back(t);
                        break;
                    }
                }

                while(stack.size() > 1 + reduceUntil &&
                      Precedence(t) <= Precedence(stack[stack.size() - 2])) {
                    if(!Reduce(error)) return false;
                }

                stack.push_back(t);
                break;

            case TokenType::UNARY_OP:
            case TokenType::OPERAND:
                stack.push_back(t);
                break;
        }
    }

    return true;
}

Expr *ExprParser::Parse(const std::string &input, std::string *error) {
    ExprParser parser;
    parser.it  = input.cbegin();
    parser.end = input.cend();
    if(!parser.Parse(error)) return NULL;

    Token r = parser.PopOperand(error);
    if(r.IsError()) return NULL;
    return r.expr;
}

} // namespace SolveSpace
