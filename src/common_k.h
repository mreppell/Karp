#ifndef KARP_COMMON_H
#define KARP_COMMON_H

#define KARP_VERSION "0.0.1"

#include <string>
#include <vector>
#include <iostream>

struct ProgramOptions {
  int threads;
  int k;
  std::string index;
  int skip; 
  int min_range;
  double em_converge;
  double illumina_version;
  std::vector<std::string> transfasta;
  std::vector<std::string> files;
  std::vector<bool> file_direction;
  bool write_index;
  bool single_end;
  bool strand_specific;
  std::string gfa; // used for inspect
  double minimum_frequency_cutoff;
  int max_em_iterations;
  //int em_restarts;
  std::vector<std::string> taxonomy_files;
  //std::string format;
  std::string out;
  bool harp_filter;
  int resc_thresh;
  bool multi_match;
  double harp_filter_thresh;
  int it_skip;
  bool collapse;
  double min_logl;
  bool fail;
  bool likeplot;
  std::string version;
  bool readinfo;

ProgramOptions() :
  min_logl(-150),
    collapse(false),
    harp_filter(false),
    threads(1),
    k(31),
    resc_thresh(500),
    skip(1),
    illumina_version(33), 
    min_range(1),
    write_index(false),
    single_end(false),
    strand_specific(false),
    minimum_frequency_cutoff(-99),
    max_em_iterations(2500)
//    em_restarts(3)
  {}
};

std::string pretty_num(size_t num);
std::string pretty_num(unsigned int num);
std::string pretty_num(int num);

void split(const std::string& str, const std::string& delimiters, std::vector<std::string>& tokens);

#endif // KARP_COMMON_H
