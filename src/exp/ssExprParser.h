#pragma once

#include <string>
#include <vector>

#include "ssExpr.h"

namespace SolveSpace
{
    //-----------------------------------------------------------------------------
// A parser; convert a string to an expression. Infix notation, with the
// usual shift/reduce approach. I had great hopes for user-entered eq
// constraints, but those don't seem very useful, so right now this is just
// to provide calculator type functionality wherever numbers are entered.
//-----------------------------------------------------------------------------


class ExprParser {
public:
    enum class TokenType {
        ERROR = 0,

        PAREN_LEFT,
        PAREN_RIGHT,
        BINARY_OP,
        UNARY_OP,
        OPERAND,

        END,
    };

    class Token {
    public:
        TokenType  type;
        Expr      *expr;

        static Token From(TokenType type = TokenType::ERROR, Expr *expr = NULL);
        static Token From(TokenType type, Expr::Op op);
        bool IsError() const { return type == TokenType::ERROR; }
    };

    std::string::const_iterator it, end;
    std::vector<Token> stack;

    char ReadChar();
    char PeekChar();

    std::string ReadWord();
    void SkipSpace();

    Token PopOperator(std::string *error);
    Token PopOperand(std::string *error);

    int Precedence(Token token);
    Token LexNumber(std::string *error);
    Token Lex(std::string *error);
    bool Reduce(std::string *error);
    bool Parse(std::string *error, size_t reduceUntil = 0);

    static Expr *Parse(const std::string &input, std::string *error);
};

} // namespace SolveSpace
