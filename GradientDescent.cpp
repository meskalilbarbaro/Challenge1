#include "GradientDescent.hpp"
#include <functional>
#include <cmath>
//Auxiliary function to calculate vector differences
std::vector<double> diff(std::vector<double> x, std::vector<double> y){
    for(size_t i=0; i<x.size(); i++){
        x[i]=x[i]-y[i];
    }
    return x;
}
//Auxiliary function to calculate vector norms
double norm(std::vector<double> x){
    double norm=0;
    for(size_t i=0; i<x.size(); i++){
        norm = norm + std::pow(x[i],2);
    }
    return std::pow(norm,0.5);
}
//Auxiliary function to calculate product by a scalar
std::vector<double> prodbyScalar(double a, std::vector<double> V){
    for(size_t i=0; i<V.size(); i++){
        V[i]=V[i]*a;
    }
    return V;
}
//Function that calculates the optimal alpha based on the Armijo rule
double Armijo(double alpha0, const double mu,const unsigned int k, const std::vector<double>& xk, FunctionWrap f, MultFunctionWrap df){
    while(f(xk) - f(diff(xk,prodbyScalar(alpha0,df(xk)))) < mu*alpha0*norm(df(xk))*norm(df(xk)) ){
        alpha0=alpha0/2;
    }
    return alpha0;
}
std::vector<double> GradientDescent(const std::vector<double>& initial_x, 
                                    const unsigned int max_it, 
                                    const double tol_df, 
                                    const double tol_x, 
                                    FunctionWrap f,
                                    MultFunctionWrap df,
                                    const double alpha0,
                                    const double mu,
                                    const int method){
    bool converged=false;
    unsigned int it =0;
    std::vector<double> xk(initial_x);


    //Construction of the Alphas vector. Each lambda function stored represents a method to evaluate optimal alphas.
    std::vector<std::function<double(const unsigned int)>> Alphas;
    Alphas.push_back([=](const unsigned int k){ return alpha0*std::exp(-mu*k);});
    Alphas.push_back([=](const unsigned int k){ return alpha0/(1+mu*k);});
    Alphas.push_back([=](const unsigned int k){ return Armijo(alpha0,mu,k,xk,f,df);});

    while(!converged){
        //Compute next xk
        std::vector<double> xk1(xk.size());
        double a=Alphas[method](it);
        std::vector<double> DF=df(xk);
        xk1=diff(xk, prodbyScalar(a,DF));
        //Check convergence
        if(it > max_it or norm(diff(xk,xk1))<tol_x or norm(df(xk)) < tol_df)
            converged=true;
        xk=xk1;
        it++;
    }

    return xk;
}
