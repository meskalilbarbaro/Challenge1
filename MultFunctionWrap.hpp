
#ifndef CHALLENGE1PACS_MULTFUNCTIONWRAP_HPP
#define CHALLENGE1PACS_MULTFUNCTIONWRAP_HPP
#include "FunctionWrap.hpp"
class MultFunctionWrap{
private:
    std::vector<FunctionWrap> function_wrap_v;
public:
    //constructor
    MultFunctionWrap(const std::vector<std::string> &expr_v, size_t size);
    //Evaluation operator
    std::vector<double> operator()(const std::vector<double> &x);
};
#endif //CHALLENGE1PACS_MULTFUNCTIONWRAP_HPP
