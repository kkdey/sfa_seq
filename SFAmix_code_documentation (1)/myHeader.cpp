#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <Eigen/Dense>
#include <limits>
#include <stdio.h>
#include <ostream>
#include <string>
#include <iostream>
#include <float.h>
#include <boost/math/special_functions/gamma.hpp>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

const double PI  =3.141592653589793238462;

void inv_psi(Eigen::MatrixXd& psi,Eigen::MatrixXd& psi_inv,int n){
    for(int i=0;i<n;i++){
        psi_inv(i,i)=double(1)/psi(i,i);
    }
}
double log_norm(double x,double e, double v){
  //const double PI  =3.141592653589793238462;
    double a = std::numeric_limits<double>::infinity();
    //cout << "a " << a << endl;
    if(v==0){
        return a;
    }else{
        return -0.5*(x-e)*(x-e)/v-0.5*log(v)-0.5*log(2*PI);
    }   
}
double log_gamma(double x, double a, double beta){
    //return a*log(beta)-log(gsl_sf_gamma(a))+(a-1)*log(x)-beta*x;
    return a*log(beta)-boost::math::lgamma(a)+(a-1)*log(x)-beta*x;
    //boost::lambda
}
void cumsum(Eigen::VectorXd& S, Eigen::VectorXd& D, int n){
  for(int i=0;i<n;i++){
    if(i=0){
      D(i)=S(i);
    }
    D(i)=D(i-1)+S(i);
  }
}

double fx(double x,double c){
    return -1*log(x)+0.5/x+c;
}

double dfx(double x){
    return -1*(double(1)/x+0.5/x/x);
}

double NR(double c){
    double x=1e-10;
    for(int i=0;i<500;i++){
        x=x-fx(x,c)/dfx(x);
    }
    return x;
}

double like(Eigen::MatrixXd& PSI,Eigen::MatrixXd& EXX,Eigen::MatrixXd& LAM,Eigen::MatrixXd& THETA,Eigen::MatrixXd& DELTA,Eigen::VectorXd& PHI,Eigen::VectorXd& TAU,Eigen::MatrixXd& Z,Eigen::MatrixXd& V,double ETA, double GAMMA,double alpha, double beta, int n,int p, int nf,double a, double b, double c, double d, double e, double f, double nu){
  double det_psi=0;
  
  for(int i=0;i<n;i++){
    det_psi = det_psi+log(PSI(i,i));
  }
  
  double like=(-1)*0.5*n*p*log(2*PI)-0.5*p*det_psi;
  
  double sum_x=0;
  for(int i=0;i<nf;i++){
    sum_x = sum_x + (-1)*0.5*EXX(i,i);
  }
  
  like = like - 0.5*nf*p*log(2*PI) + sum_x;
  
  for(int i=0;i<nf;i++){
    like=like + (Z(0,i)+alpha)*V(0,i) + (Z(1,i)+beta)*V(1,i);
  }
  
  for(int i=0;i<n;i++){
    for(int j=0;j<nf;j++){
      if(THETA(i,j)!=0){
        like=like+Z(0,j)*log_norm(LAM(i,j),0,THETA(i,j));
        //if(DELTA(i,j)!=0){
          like=like+Z(0,j)*log_gamma(THETA(i,j),a,DELTA(i,j));
          //}
      }
      like=like+Z(0,j)*log_gamma(DELTA(i,j),b,PHI(j));
      like=like+Z(1,j)*log_norm(LAM(i,j),0,PHI(j));
    }
  }
 
 for(int i=0;i<nf;i++){
   like=like+log_gamma(PHI(i),c,TAU(i));
   like=like+log_gamma(TAU(i),d,ETA);
 }
 
 like=like+log_gamma(ETA,e,GAMMA);
 like=like+log_gamma(GAMMA,f,nu);
  
  return like;
 
}
/*
#include <iostream>     // cout
#include <math.h>       // acos
#include <float.h>      // DBL_MAX
#include <limits>       // numeric_limits
*/
template<typename T>
bool is_infinite( const T &value )
{
    // Since we're a template, it's wise to use std::numeric_limits<T>
    //
    // Note: std::numeric_limits<T>::min() behaves like DBL_MIN, and is the smallest absolute value possible.
    //
 
    T max_value = std::numeric_limits<T>::max();
    T min_value = - max_value;
 
    return ! ( min_value <= value && value <= max_value );
}
 
template<typename T>
bool is_nan( const T &value )
{
    // True if NAN
    return value != value;
}
 
template<typename T>
bool is_valid( const T &value )
{
    return ! is_infinite(value) && ! is_nan(value);
}

/*
int main()
{
    using std::cout;
 
    double a, b, c, d, e;
 
    a = 1.0;
    b = 0.0;
    c = a / c;          // divide by zero
    d = acos(-1.001);   // domain for acos is [-1, 1], anything else is #IND or inf
    e = b / b;          // zero / zero
 
    cout << "Value of a: " << a << " " << is_valid(a) << " " << (is_nan(a) ? " nan " : "") << (is_infinite(a) ? " infinite " : "") << "n";
    cout << "Value of b: " << b << " " << is_valid(b) << " " << (is_nan(b) ? " nan " : "") << (is_infinite(b) ? " infinite " : "") << "n";
    cout << "Value of c: " << c << " " << is_valid(c) << " " << (is_nan(c) ? " nan " : "") << (is_infinite(c) ? " infinite " : "") << "n";
    cout << "Value of d: " << d << " " << is_valid(d) << " " << (is_nan(d) ? " nan " : "") << (is_infinite(d) ? " infinite " : "") << "n";
    cout << "Value of e: " << e << " " << is_valid(e) << " " << (is_nan(e) ? " nan " : "") << (is_infinite(e) ? " infinite " : "") << "n";
 
    return 0;
}
*/
