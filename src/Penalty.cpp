#include "RcppArmadillo.h"
#include "Penalty.h"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]


// =============================================================================
// Proximal operators

arma::rowvec proximal_L1(
    const arma::rowvec &b,
    const arma::rowvec &m
){
  uint nt = b.n_elem;
  arma::rowvec s(nt);
  s = arma::abs(b) - m;
  s = arma::sign(b) % arma::clamp(s, 0., arma::datum::inf);
  return s;
}

arma::rowvec proximal_L2(
    const arma::rowvec &s,
    const double m
){
  double sn = arma::norm(s);
  if(sn == 0.) sn = 1.; // to avoid dividing by 0, the value doesn't matter since s=0 in that case
  return fmax(1. - m/sn, 0.) * s;
}

arma::rowvec proximal_L1L2_row(
    const arma::rowvec &b,
    const arma::rowvec &m1,
    const double m2
){
  return proximal_L2(proximal_L1(b, m1), m2);
}

arma::mat proximal_L1L2(
    const arma::mat &B,
    const arma::mat &M1,
    const arma::colvec &m2
){
  arma::mat out = B * 0.;
  for(uint j=0; j<B.n_rows; j++){
    out.row(j) = proximal_L1L2_row(B.row(j), M1.row(j), m2(j));
  }
  return out;
}

// =============================================================================
// Interface + Factory
Penalty* Penalty::Create(
    std::string name,
    double lambda,
    double alpha,
    double power,
    double a,
    bool penalize_intercept
){
  if(name=="lasso") return new AdaptiveSparseGroupLasso(lambda, 1., 0., penalize_intercept);
  if(name=="group_lasso") return new AdaptiveSparseGroupLasso(lambda, 0., 0., penalize_intercept);
  if(name=="sparse_group_lasso") return new AdaptiveSparseGroupLasso(lambda, alpha, 0., penalize_intercept);
  if(name=="adaptive_lasso") return new AdaptiveSparseGroupLasso(lambda, 1., power, penalize_intercept);
  if(name=="adaptive_group_lasso") return new AdaptiveSparseGroupLasso(lambda, 0., power, penalize_intercept);
  if(name=="adaptive_sparse_group_lasso") return new AdaptiveSparseGroupLasso(lambda, alpha, power, penalize_intercept);

  if(name=="bridge") return new AdaptiveSparseGroupLasso(lambda, 1., 1.-power, penalize_intercept);
  if(name=="group_bridge") return new AdaptiveSparseGroupLasso(lambda, 0., 1.-power, penalize_intercept);
  if(name=="sparse_group_bridge") return new AdaptiveSparseGroupLasso(lambda, alpha, 1.-power, penalize_intercept);

  if(name=="scad") return new SparseGroupSCAD(lambda, 1., a, penalize_intercept);
  if(name=="group_scad") return new SparseGroupSCAD(lambda, 0., a, penalize_intercept);
  if(name=="sparse_group_scad") return new SparseGroupSCAD(lambda, alpha, a, penalize_intercept);

  if(name=="mcp") return new SparseGroupMCP(lambda, 1., a, penalize_intercept);
  if(name=="group_mcp") return new SparseGroupMCP(lambda, 0., a, penalize_intercept);
  if(name=="sparse_group_mcp") return new SparseGroupMCP(lambda, alpha, a, penalize_intercept);
  Rcpp::stop("unrecognized penalty name: " + name);
}

double Penalty::eval(const arma::mat &B){
  double l1term = arma::accu(this->W1 % arma::abs(B));
  double l2term = arma::accu(this->w2 % arma::sqrt(arma::sum(arma::pow(B, 2), 1)));
  return this->alpha * l1term + (1 - this->alpha) * l2term;
}

void Penalty::update_B0(const arma::mat &B){
  this->AbsB0 = arma::abs(B);
  this->normB0 = arma::sqrt(arma::sum(arma::pow(B, 2), 1));
}

void Penalty::update_weights(){
  this->W1 = this->derivative(this->AbsB0, 1.);
  this->w2 = this->derivative(this->normB0, sqrt(this->AbsB0.n_cols));
  if(!this->penalize_intercept){
    this->W1.row(0).zeros();
    this->w2(0) = 0.;
  }
}

void Penalty::unit_weights(const arma::mat &B){
  // without lambda
  this->W1 = arma::mat(B.n_rows, B.n_cols, arma::fill::ones);
  this->w2 = arma::colvec(B.n_rows, arma::fill::ones) * sqrt(B.n_cols);
  if(!this->penalize_intercept){
    this->W1.row(0).zeros();
    this->w2(0) = 0.;
  }
}

void Penalty::large_weights(const arma::mat &B){
  this->unit_weights(B);
  this->W1 *= 1e10;
  this->w2 *= 1e10;
}

arma::mat Penalty::proximal(const arma::mat &B, double stepsize){
  arma::mat M1 = this->W1 * this->alpha * stepsize;
  arma::colvec m2 = this->w2 * (1 - this->alpha) * stepsize;
  return proximal_L1L2(B, M1, m2);
}

double Penalty::lambda_max(arma::mat B, const arma::mat &gB, const double stepsize){
  // TODO: check if W1, w2 are fine here or it should be AbsB0, normB0
  // NB all penalties have the same lambda_max (I think)
  // basically since they have linear penalty at zero, so SGL is fine
  // eq. derivative at zero is always lambda
  // gB.print("gB for lambda max");
  // this->W1.print("W1 for lambda max");
  // this->w2.print("w2 for lambda max");
  double lambda_max = 0.;
  double num, denum, gnorm, wnorm;
  for(uint j=0; j<gB.n_rows; j++){
    if(j==0 and !this->penalize_intercept) continue;
    // L1 bound
    if(this->alpha > 0.){
      for(uint t=0; t<gB.n_cols; t++){
        if(this->W1(j, t) > 0.){
          num = fabs(gB(j, t));
          denum = this->alpha * this->W1(j, t) + (1-this->alpha) * this->w2(j);
          lambda_max = fmax(lambda_max, num / denum);
        }
      }
    }
    // L2 bound
    if(this->alpha < 1. and this->w2(j) > 0.){
      gnorm = arma::norm(gB.row(j), "fro");
      wnorm = arma::norm(this->W1.row(j), "fro");
      lambda_max = fmax(lambda_max, gnorm / (this->alpha * wnorm + (1-this->alpha) * this->w2(j)));
    }
  }
  // backtrack
  arma::mat tmpB = B;
  // this->W1.print("W1");
  // tmpB.print("tmpB");
  double factor = 1.;
  this->lambda = lambda_max;
  this->update_weights();
  while(arma::accu(this->W1 % arma::abs(tmpB)) < 1e-10){
    factor /= 1.1;
    tmpB = B - gB * stepsize;
    tmpB = proximal_L1L2(
      tmpB,
      this->W1 * this->alpha * factor * stepsize,
      this->w2 * (1 - this->alpha) * factor * stepsize
    );
    // Rcpp::Rcout << "lambda_max : " << lambda_max << std::endl;
    // tmpB.print("tmpB");
  }
  factor *= 1.1; // when we escape the while loop, we have gone too far
  // Rcpp::Rcout << "lambda_max : " << lambda_max << std::endl;
  // Rcpp::Rcout << "factor : " << factor << std::endl;
  return lambda_max*factor;
}



// =============================================================================
// Adaptive SGL

AdaptiveSparseGroupLasso::AdaptiveSparseGroupLasso(
  double lambda,
  double alpha,
  double power,
  bool penalize_intercept
){
  this->lambda = lambda;
  this->alpha = alpha;
  this->power = power;
  this->penalize_intercept = penalize_intercept;
  if(power==0.) this->compute_penalty_weights = false;
}

void AdaptiveSparseGroupLasso::add_to_results(Rcpp::List& results){
  results["penalty.name"] = "adaptive_sparse_group_lasso";
  results["penalty.alpha"] = this->alpha;
  results["penalty.lambda"] = this->lambda;
  results["penalty.adaptive"] = this->power;
  results["penalty.penalize_intercept"] = this->penalize_intercept;
}

arma::mat AdaptiveSparseGroupLasso::derivative(const arma::mat &W, double scaling){
  arma::mat Bpow = arma::pow(W, -this->power);
  Bpow.elem(arma::find(W < 1e-10)).fill(1e10);
  arma::mat out = scaling * this->lambda * Bpow;
  return out;
}


// =============================================================================
// SCAD

SparseGroupSCAD::SparseGroupSCAD(
  double lambda,
  double alpha,
  double a,
  bool penalize_intercept
){
  this->lambda = lambda;
  this->alpha = alpha;
  this->a = a;
  this->penalize_intercept = penalize_intercept;
}

void SparseGroupSCAD::add_to_results(Rcpp::List& results){
  results["penalty.name"] = "sparse_group_scad";
  results["penalty.alpha"] = this->alpha;
  results["penalty.lambda"] = this->lambda;
  results["penalty.a"] = this->a;
  results["penalty.penalize_intercept"] = this->penalize_intercept;
}

arma::mat SparseGroupSCAD::derivative(const arma::mat &W, double scaling){
  double lam = this->lambda * scaling;
  arma::mat cond = arma::conv_to<arma::mat>::from(W < lam); // convert boolean to double
  arma::mat out = lam * cond;
  out += (1. - cond) % arma::clamp(this->a * lam - W, 0., arma::datum::inf) / (this->a - 1.);
  // arma::mat out = scaling * this->lambda * arma::pow(W, -0.5);
  return out;
}


// =============================================================================
// MCP

SparseGroupMCP::SparseGroupMCP(
  double lambda,
  double alpha,
  double a,
  bool penalize_intercept
){
  this->lambda = lambda;
  this->alpha = alpha;
  this->a = a;
  this->penalize_intercept = penalize_intercept;
}

void SparseGroupMCP::add_to_results(Rcpp::List& results){
  results["penalty.name"] = "sparse_group_mcp";
  results["penalty.alpha"] = this->alpha;
  results["penalty.lambda"] = this->lambda;
  results["penalty.a"] = this->a;
  results["penalty.penalize_intercept"] = this->penalize_intercept;
}

arma::mat SparseGroupMCP::derivative(const arma::mat &W, double scaling){
  double lam = this->lambda * scaling;
  arma::mat cond = arma::conv_to<arma::mat>::from(W <= lam * this->a); // convert boolean to double
  arma::mat out = (lam - W / this->a) % cond;
  return out;
}
