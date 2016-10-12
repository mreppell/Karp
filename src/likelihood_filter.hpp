#ifndef KARP_LIKEFILTER
#define KARP_LIKEFILTER

#include<zlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <karpeigen/Eigen/Core>
#include <karpeigen/Eigen/Dense>

#include "common_k.h"

using Eigen::ArrayXXi;

struct LKFilter {

  double mean;
  double sd;
  bool use;
  int numreads;
  bool paired;
  const double& illumina_version;
  double threshold;
  unsigned int number_removed;

  LKFilter(const ProgramOptions& opt) : use(false), mean(0.0), sd(0.0), illumina_version(opt.illumina_version), paired(opt.single_end), numreads(0), threshold(opt.harp_filter_thresh), number_removed(0) {}

  void GetParameterValues(const ProgramOptions& opt);
  std::vector<double> E_logl(ArrayXXi& hf_qualityscores);
  std::vector<double> E_logl_per_site(size_t position,ArrayXXi& hf_qualityscores);
  std::vector<double> E_logl_per_site_conditional(size_t q);
  void lk_processBuffer(ArrayXXi& hf_qualityscores,std::vector<std::pair<const char*, int>>& quals);
  void updateFilter(ArrayXXi& hf_qualityscores,const char* quals1,const char* quals2,bool paired, const double& illumina_version);
  void justGetReadNum(const ProgramOptions& opt);

};

#endif // KARP_LIKEFILTER
