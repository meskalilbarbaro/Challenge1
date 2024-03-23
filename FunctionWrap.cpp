#include "FunctionWrap.hpp"
FunctionWrap::FunctionWrap(const std::string & expr, size_t size){
    inputs.resize(size);
    for(size_t i=0;i<size;i++){
        function_parser.DefineVar("x"+std::to_string(i+1), &inputs[i]);
    }
    function_parser.SetExpr(expr);
}
double FunctionWrap::operator()(const std::vector<double> &x){
    inputs=x;
    return function_parser.Eval();
}
FunctionWrap::FunctionWrap(const FunctionWrap & f1): function_parser(f1.function_parser){
    inputs.resize(f1.inputs.size());
    for(size_t i=0;i<inputs.size();i++){
        function_parser.DefineVar("x"+std::to_string(i+1), &inputs[i]);
    }
}