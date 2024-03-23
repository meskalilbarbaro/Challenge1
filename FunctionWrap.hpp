#ifndef CHALLENGE1PACS_FUNCTIONWRAP_HPP
#define CHALLENGE1PACS_FUNCTIONWRAP_HPP
#include <vector>
#include <muParser.h>
#include <memory>
#include <string>
class FunctionWrap{
private:
    std::vector<double> inputs;
    mu::Parser function_parser;
public:
    FunctionWrap(const std::string & expr, size_t size);
    FunctionWrap(const FunctionWrap & f1);
    double operator()(const std::vector<double> &x);
};
#endif //CHALLENGE1PACS_FUNCTIONWRAP_HPP
