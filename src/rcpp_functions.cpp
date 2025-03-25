#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Erlang_dist(NumericVector  x, int alpha, double beta){
  return exp(-x/beta-lgamma(alpha)+(alpha-1)*log(x))/pow(beta,alpha);
}

// [[Rcpp::export]]
NumericVector Erlang_mix2(NumericVector x, NumericVector support, NumericVector weight, double theta){
  NumericVector tmp(x.length());
  int p=support.length();
  for (int i=0;i<p;i++){
    tmp +=weight(i)*Erlang_dist(x,support(i),theta);
  }
  return tmp;
}

// [[Rcpp::export]]
NumericVector Erlang_mix(NumericVector x, NumericVector support, NumericVector weight, double theta){
  NumericVector tmp(x.length());
  int p=support.length();
  for (int i=0;i<p;i++){
    tmp +=weight(i)*exp(-x/theta-lgamma(support(i))+(support(i)-1)*log(x))/pow(theta,support(i));
  }
  return tmp;
}


// [[Rcpp::export]]
NumericVector mix_ratio(NumericVector x, int phi, NumericVector support, NumericVector weight, double theta){
  NumericVector tmp(x.length());
  NumericVector numer=exp(-lgamma(phi)+(phi-1)*log(x))/pow(theta,phi);
  int p=support.length();
  for (int i=0;i<p;i++){
    tmp +=weight(i)*exp(-lgamma(support(i))+(support(i)-1)*log(x))/pow(theta,support(i));
  }
  return numer/tmp;
}


// [[Rcpp::export]]
double Erlang_gradient(NumericVector x, int phi, NumericVector support, NumericVector weight, double theta){
  return sum(mix_ratio(x,phi,support,weight,theta))-x.length();
}

// [[Rcpp::export]]
double Erlang_dgradient(NumericVector x, int phi, NumericVector support, NumericVector weight, double theta){
  return sum((log(x)-log(theta)-R::digamma(phi))*mix_ratio(x,phi,support,weight,theta));
}

// [[Rcpp::export]]
double Erlang_d2gradient(NumericVector x, int phi, NumericVector support, NumericVector weight, double theta){
  return sum((pow(log(x)-log(theta)-R::digamma(phi),2)-R::trigamma(phi))*mix_ratio(x,phi,support,weight,theta));
}

// [[Rcpp::export]]
NumericMatrix bbb(NumericVector x, NumericVector support, NumericVector weight, double theta){
  NumericMatrix tmp(x.length(),support.length());
  for (int j=0;j<support.length();j++){
    tmp(_,j)=Erlang_dist(x,support(j),theta);
  }
  return tmp;
}

// [[Rcpp::export]]
NumericMatrix ccc(NumericMatrix tmp, NumericVector weight){
  NumericMatrix tmp11(tmp.nrow(),tmp.ncol());
  for (int j=0;j<weight.length();j++){
    tmp11(_,j)=weight(j)*tmp(_,j);
  }
  return tmp11;
}


// [[Rcpp::export]]
NumericMatrix compute_S(NumericVector x, NumericVector support, NumericVector weight, double theta){
  int p=support.length();
  int n=x.length();
  NumericMatrix tmp(n,p);
  NumericMatrix val(n,p);
  NumericMatrix tmp11(n,p);
  for (int j=0;j<p;j++){
    tmp(_,j)=Erlang_dist(x,support(j),theta);
    tmp11(_,j)=weight(j)*tmp(_,j);
  }
  for (int i=0;i<n;i++){
    val(i,_)=tmp(i,_)/Rcpp::sum(tmp11(i,_));
  }
  return val;
}



// [[Rcpp::export]]
NumericVector check_descent(NumericVector x, NumericVector support, NumericVector weight, NumericVector new_weight, double theta,double lik){
  int count=0;
  while (lik>sum(log(Erlang_mix(x,support,new_weight,theta)))){
    count++;
    new_weight=weight+pow(0.5,count)*(new_weight-weight);
    new_weight=new_weight/Rcpp::sum(new_weight);
    if (count==15){
      new_weight=weight;
      break;
    }
  }
  return new_weight;
}










