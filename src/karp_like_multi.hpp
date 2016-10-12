#ifndef KARP_LIKE_MULTI_H
#define KARP_LIKE_MULTI_H

#include<string>
#include<vector>
#include<unordered_map>
#include<karpeigen/Eigen/Core>
#include<karpeigen/Eigen/Dense>
#include <zlib.h>
#include <thread>
#include <sstream>


#include "fastaIndex.h"
#include "common_k.h"

using Eigen::ArrayXXd;
using Eigen::ArrayXXi;
using Eigen::ArrayXd;
using Eigen::ArrayXi;

struct thread_output {
  std::vector<std::unordered_map<double,std::vector<unsigned int> > > clikes;
  //  ArrayXXi c_hf_qualityscores;
  std::unordered_map<unsigned int,bool > c_z_support;
  unsigned int numreads;
  unsigned int numresc;
  unsigned int numstrict;
  std::string mapfail;
  std::unordered_map<int, std::string> read_order;
};

class LikelihoodProcessor {
  
public:
  LikelihoodProcessor(std::unordered_map<unsigned int,double>& current_freqs,int threads,double& convergence_threshold,double& num_single,std::unordered_map<unsigned int, double>& prior_counts,std::string& out_base,std::unordered_map<unsigned int,bool >& zero_support,double& min_logl) : threads(threads), convergence_threshold(convergence_threshold), out_base(out_base), num_single(num_single), prior_counts(prior_counts),zero_support(zero_support),min_logl(min_logl),final_logl(0.0) {
 
    std::unordered_map<unsigned int,double> pre_abundance_counts;
    for (auto get=current_freqs.begin();get!=current_freqs.end();++get) {
      std::pair<unsigned int, double> abent(get->first,0.0);
      pre_abundance_counts.insert(abent);
    }
    for (int i=0;i<threads;++i) {
      abundance_counts.push_back(pre_abundance_counts);
      read_counts.push_back(0);
      marginal_likes.push_back(0);
    }
  }
  
  double final_logl;
  double& min_logl;
  std::unordered_map<unsigned int,bool >& zero_support;
  std::string& out_base;
  double& convergence_threshold;
  int threads;
  std::vector<std::unordered_map<unsigned int,double> > abundance_counts;
  std::vector<double> read_counts;
  std::vector<double> marginal_likes;
  std::unordered_map<unsigned int,double> freqs;
  double& num_single;
  std::unordered_map<unsigned int, double>& prior_counts;
  

  void ProcessLikes(std::unordered_map<unsigned int,double>& current_freqs,double& frequency_cutoff);

};

struct SquaremControl {

  int K;
  int method;
  double mstep;
  int maxiter;
  bool square;
  bool trace;
  double stepmin0;
  double stepmax0;
  double kr;
  double objfninc;
  double tol;

  SquaremControl(ProgramOptions& opt) : K(1), method(3), mstep(4), maxiter(opt.max_em_iterations), square(true), trace(false), stepmin0(1), stepmax0(1), kr(1), objfninc(1), tol(opt.em_converge) {}

};

struct SquaremOutput{
  std::unordered_map<unsigned int, double> freqs;
  double final_logl;
  int iter;
  //uint64_t final_single_fail;
  int pfevals;
  bool convergence;
};

struct FPTOut {
  std::unordered_map<unsigned int, double> freqs;
  double final_logl;
  bool at_min;
  //uint64_t final_single_fail;
};

void GetAbundance(LikelihoodProcessor* LP,int thread_num,std::unordered_map<unsigned int,double> current_freqs,double base_marginal,std::string out_base,double min_logl);
//void GetAbundance(int thread_num,std::unordered_map<unsigned int,double> current_freqs,double base_marginal,std::string out_base,double min_logl);
void cutoff_low_frequencies(std::unordered_map<unsigned int,double>& f_old,std::unordered_map<unsigned int,double>& f_new,double& min_freq);


struct HLK {

  std::unordered_map<unsigned int,bool > zero_support;
  double min_logl;
  std::unordered_map<unsigned int,double> freqs;
  std::unordered_map<unsigned int, double> prior_counts;
  double remaining_freq;
  double convergence_threshold;
  double nummapped;
  double num_single;
  double base_marginal;
 
  
  double minimum_frequency_cutoff;
  double current_frequency_cutoff;
  double frequency_cutoff_iteration;
  std::stringstream& logfile_track;
  const double& illumina_version;
  int& threads;
  std::string& out_base;
  double num_single_fail;
  double max_em_iterations;
  SquaremControl EMcontrol;
 
  HLK(fastaIndex& findex,ProgramOptions& opt,std::stringstream& logfile_track) : convergence_threshold(opt.em_converge), EMcontrol(opt), min_logl(opt.min_logl), minimum_frequency_cutoff(opt.minimum_frequency_cutoff), out_base(opt.out), logfile_track(logfile_track), illumina_version(opt.illumina_version), threads(opt.threads), num_single(0), max_em_iterations(opt.max_em_iterations), nummapped(0), base_marginal(0), num_single_fail(0) {}
  
  bool StreamReads(std::unordered_map<unsigned int,double>& current_freqs,double& iteration,double& current_frequency_cutoff,double& frequency_cutoff_iteration);
  //std::unordered_map<unsigned int,double> GetPriors(std::stringstream& logfile_track);
void GetPriors(std::stringstream& logfile_track);
    
  void estimate_haplotype_frequencies();
  double initializeFreqs(std::unordered_map<unsigned int,double>& freqs,std::unordered_map<unsigned int,bool>& zero_support,std::unordered_map<unsigned int,double>& prior_freqs,double& nummapped);  
  //bool UpdateFrequencies(std::unordered_map<unsigned int,double>& current_freqs, std::vector< std::unordered_map<unsigned int,double > >& elks,double& forgetfactor,double& current_frequency_cutoff);
  SquaremOutput SquareEM(std::unordered_map<unsigned int, double>& current_freqs);
  FPTOut FPTfn(std::unordered_map<unsigned int,double>& current_freqs,bool get_to_min);

};

std::vector<int> interpretCigar(const char& cigar_string,double& deletions);
double getError(const char& qual,double mod);
double calcLikelihood(std::string& cigar_string, const char* quals,const double& illumina_version);
//double getMod(double& illumina_version);

#endif // KARP_LIKE_MULTI_H
