#include <iostream>
#include <vector>
#include <muParser.h>
#include "FunctionWrap.hpp"
#include "MultFunctionWrap.hpp"
#include <fstream>
#include "json.hpp"
#include "GradientDescent.hpp"
using json = nlohmann::json;
int main(){
    std::ifstream fdata("data.json");
    json data = json::parse(fdata);
    //standard input parameters
    const unsigned int max_it = data.value("max_it",100);
    const double tol_df = data.value("tol_df", 1e-7);
    const double tol_x = data.value("tol_x", 1e-7);
    const double alpha0 = data.value("alpha0", 0.1);
    const double decayParameter = data.value("decayParameter", 0.1);
    const int learningMethod = data.value("learningMethod", 0);
    //mandatory inputs
    std::string expr = data["fun"];
    std::vector<std::string> expr_v = data["dfun"];
    std::vector<double> initial_x = data["initial"];
    //
    FunctionWrap f(expr,initial_x.size());
    MultFunctionWrap df(expr_v,initial_x.size());
    ///////////////////////////////////////////////////////////////////
    std::vector<double> x_min = GradientDescent(initial_x, max_it, tol_df, tol_x, f, df, alpha0, decayParameter, learningMethod);

    std::cout << x_min[0] << " " << x_min[1] << std::endl;



    return 0;
}
