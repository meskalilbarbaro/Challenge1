#include "FunctionWrap.hpp"
#include "MultFunctionWrap.hpp"
#include <vector>

std::vector<double> GradientDescent(const std::vector<double>& initial_x, 
                                    const unsigned int max_it, 
                                    const double tol_df, 
                                    const double tol_x, 
                                    FunctionWrap f,
                                    MultFunctionWrap df,
                                    const double alpha0,
                                    const double mu,
                                    const int method);

