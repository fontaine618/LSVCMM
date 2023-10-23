#include "RcppArmadillo.h"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Logger_hpp
#define Logger_hpp

class FunctionProfile {
public:
  uint calls;
  double time;
  FunctionProfile() : calls(0), time(0) {}
  void add_call(const std::chrono::duration<double> diff){
    this->calls++;
    this->time += std::chrono::duration <double, std::milli> (diff).count();
  }
  Rcpp::List to_rcpp(){
    return Rcpp::List::create(
      Rcpp::Named("calls") = Rcpp::wrap(calls),
      Rcpp::Named("time") = Rcpp::wrap(time)
    );
  }
};

class Profiler {
  std::map <std::string, FunctionProfile> profiler;
public:
  void add_call(std::string name, const std::chrono::duration<double> diff){
    if (profiler.find(name) == profiler.end()){
      profiler[name] = FunctionProfile();
    }
    profiler[name].add_call(diff);
  }
  Rcpp::List to_rcpp(){
    Rcpp::List result;
    for (auto it = profiler.begin(); it != profiler.end(); ++it){
      result[it->first] = it->second.to_rcpp();
    }
    return result;
  }
};

class Logger {

  // per-round data
  std::vector<uint> round;
  std::vector<uint> mean_iterations;
  std::vector<uint> variance_iterations;
  std::vector<double> log_likelihood;

  // per mean iteration data
  std::vector<uint> mean_round;
  std::vector<uint> mean_iteration_within_round;
  std::vector<double> mean_objective;

  // per variance iteration data
  std::vector<uint> variance_round;
  std::vector<uint> variance_iteration_within_round;
  std::vector<double> variance_objective;

public:

  Profiler profiler;

  void add_round_results(
      uint round,
      uint mean_iterations,
      uint variance_iterations,
      double log_likelihood
  ){
    this->round.push_back(round);
    this->mean_iterations.push_back(mean_iterations);
    this->variance_iterations.push_back(variance_iterations);
    this->log_likelihood.push_back(log_likelihood);
  }

  void add_mean_iteration_results(
      uint round,
      uint mean_iteration_within_round,
      double mean_objective
  ){
    this->mean_round.push_back(round);
    this->mean_iteration_within_round.push_back(mean_iteration_within_round);
    this->mean_objective.push_back(mean_objective);
  }

  void add_variance_iteration_results(
      uint round,
      uint variance_iteration_within_round,
      double variance_objective
  ){
    this->variance_round.push_back(round);
    this->variance_iteration_within_round.push_back(variance_iteration_within_round);
    this->variance_objective.push_back(variance_objective);
  }

  void reset(){
    round.clear();
    mean_iterations.clear();
    variance_iterations.clear();
    log_likelihood.clear();
    mean_round.clear();
    mean_iteration_within_round.clear();
    mean_objective.clear();
    variance_round.clear();
    variance_iteration_within_round.clear();
    variance_objective.clear();
  }

  Rcpp::List to_rcpp(){
    return Rcpp::List::create(
      Rcpp::Named("round") = Rcpp::wrap(round),
      Rcpp::Named("mean_iterations") = Rcpp::wrap(mean_iterations),
      Rcpp::Named("variance_iterations") = Rcpp::wrap(variance_iterations),
      Rcpp::Named("log_likelihood") = Rcpp::wrap(log_likelihood),
      Rcpp::Named("mean_round") = Rcpp::wrap(mean_round),
      Rcpp::Named("mean_iteration_within_round") = Rcpp::wrap(mean_iteration_within_round),
      Rcpp::Named("mean_objective") = Rcpp::wrap(mean_objective),
      Rcpp::Named("variance_round") = Rcpp::wrap(variance_round),
      Rcpp::Named("variance_iteration_within_round") = Rcpp::wrap(variance_iteration_within_round),
      Rcpp::Named("variance_objective") = Rcpp::wrap(variance_objective),
      Rcpp::Named("profiler") = profiler.to_rcpp()
    );
  }
};


#endif
