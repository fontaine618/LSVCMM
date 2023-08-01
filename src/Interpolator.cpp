#include "RcppArmadillo.h"
#include "Interpolator.h"

//[[Rcpp::depends(RcppArmadillo)]]

Interpolator::Interpolator(const arma::rowvec &time){
  this->time = time;
}

arma::mat Interpolator::interpolator_matrix(const arma::colvec &new_time){
  // expect newtime n x 1
  // return n x nt
  uint n = new_time.n_rows;
  uint nt = this->time.n_cols;
  arma::mat out = arma::zeros(n, nt);
  for(uint i=0; i<n; i++){
    double ti = new_time(i);
    uint which = 0;
    for(uint j=0; j<nt; j++){
      if(this->time(j) >= ti){
        which = j;
        break;
      }
      which = nt;
    }
    // if which = nt, then ti is larger than everything
    if(which >= nt){
      out(i, nt-1) = 1.;
      continue;
    }
    // if which = 0, then ti is smaller than everything
    if(which <= 0){
      out(i, 0) = 1.;
      continue;
    }
    // otherwise: linear interpolation
    double t0 = this->time(which-1);
    double t1 = this->time(which);
    double dt = (t1 - t0);
    double t = (ti - t0) / dt;
    out(i, which-1) = 1 - t;
    out(i, which) = t;
    // cubic interpolation with finite diffs
    // double h00 = (1.+2.*t)*(1.-t)*(1.-t);
    // double h10 = t*(1.-t)*(1.-t);
    // double h01 = t*t*(3.-2.*t);
    // double h11 = t*t*(t-1.);
    // // if which == 1, we cannot do finite diff on the left
    // if(which==1){
    //   double t2 = this->time(which+1);
    //   double dt2 = (t2 - t1);
    //   // out(i, which-1) = h00 - 0.5*h11;
    //   // out(i, which) = h01 + 0.5*dt*h11*(1./dt-1./dt2);
    //   out(i, which-1) = h00 - 0.5*h11 - h10;
    //   out(i, which) = h01 + 0.5*dt*h11*(1./dt-1./dt2) + h10;
    //   out(i, which+1) = 0.5*dt*h11/dt2;
    //   continue;
    // }
    // // if which == nt-2, we cannot do finite diff on the right
    // if(which==nt-1){
    //   double tm = this->time(which-2);
    //   double dt0 = (t0 - tm);
    //   out(i, which-2) = -0.5*dt*h10/dt0;
    //   // out(i, which-1) = h00 + 0.5*dt*h10*(1./dt0-1./dt);
    //   // out(i, which) = h01 + 0.5*h10;
    //   out(i, which-1) = h00 - h11 + 0.5*dt*h10*(1./dt0-1./dt);
    //   out(i, which) = h01 + h11 + 0.5*h10;
    //   continue;
    // }
    // // cubic interpolation with finite diff tangent
    // double t2 = this->time(which+1);
    // double tm = this->time(which-2);
    // double dt2 = (t2 - t1);
    // double dt0 = (t0 - tm);
    // out(i, which-2) = -0.5*dt*h10/dt0;
    // out(i, which-1) = h00 - 0.5*h11 + 0.5*dt*h10*(1./dt0-1./dt);
    // out(i, which) = h01 + 0.5*dt*h11*(1./dt-1./dt2) + 0.5*h10;
    // out(i, which+1) = 0.5*dt*h11/dt2;
  }
  return out;
}

std::vector<arma::mat> Interpolator::interpolator_matrix(const std::vector<arma::colvec> &new_time){
  std::vector<arma::mat> out(new_time.size());
  for(uint i=0; i<new_time.size(); i++) out[i] = this->interpolator_matrix(new_time[i]);
  return out;
}
