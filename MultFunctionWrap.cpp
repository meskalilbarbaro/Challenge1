#include "MultFunctionWrap.hpp"
MultFunctionWrap::MultFunctionWrap(const std::vector<std::string> &expr_v, size_t size){
    for(const std::string & expr : expr_v){

        function_wrap_v.emplace_back(expr,size);
    }
}
std::vector<double> MultFunctionWrap::operator()(const std::vector<double> &x){
    std::vector<double> result;
    for(FunctionWrap f : function_wrap_v){


        result.emplace_back(f(x));
    }
    return result;
}