#ifndef KARP_PROCESSREADS_H
#define KARP_PROCESSREADS_H

#include <zlib.h>
#include "kseq.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>

#include "MinCollector.h"
#include "ssw_cpp.h"
#include "Kmer.hpp"
#include "common_k.h"
#include "karp_like_multi.hpp"
#include "likelihood_filter.hpp"
#include "taxonomy.hpp"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

uint32_t ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc, fastaIndex& findex,std::stringstream& logfile_track,LKFilter& like_filter,EarlyTaxonomy& earlytax,HLK& likelihoods);

class SequenceReader {
public:

  SequenceReader(const ProgramOptions& opt) :
  fp1(0),fp2(0),seq1(0),seq2(0),
  l1(0),l2(0),nl1(0),nl2(0),
  paired(opt.single_end), files(opt.files),
    current_file(0), state(false) {}
  bool empty();
  ~SequenceReader();

  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals);

private:
  gzFile fp1 = 0, fp2 = 0;
  kseq_t *seq1 = 0, *seq2 = 0;
  int l1,l2,nl1,nl2;
  bool paired;
  const std::vector<std::string>& files;
  int current_file;
  bool state; // is the file open
};

class MasterProcessor {
public:
  MasterProcessor (KmerIndex &index, const ProgramOptions& opt, MinCollector &tc, fastaIndex& findex, LKFilter& like_filter,std::stringstream& logfile_track, EarlyTaxonomy& earlytax)
    : tc(tc), index(index), opt(opt), SR(opt), numreads(opt.threads,0), numresc(opt.threads,0), numstrict(opt.threads,0), zsupport(opt.threads), findex(findex), logfile_track(logfile_track), like_filter(like_filter), earlytax(earlytax), collapse(opt.collapse) {} 
 
  std::mutex reader_lock;
  std::mutex writer_lock;

  LKFilter& like_filter;
  EarlyTaxonomy& earlytax;
  std::stringstream& logfile_track;
  fastaIndex& findex;
  SequenceReader SR;
  MinCollector& tc;
  KmerIndex& index;
  const ProgramOptions& opt;
  std::vector<uint32_t> numreads;
  std::vector<uint32_t> numstrict;
  std::vector<uint32_t> numresc;
  std::vector<std::unordered_map<unsigned int, bool> > zsupport;
  bool collapse;

  void processReads();

  void update(const std::vector<int>& c, const std::vector<std::vector<int>>& newEcs, int n,thread_output thread_out);
};

struct rt_entry {
  std::unordered_map<double, int> counts;
  std::unordered_map<std::string, double> seen_likes;
  std::vector<unsigned int> single_starts;
  std::vector<unsigned int> paired_starts;

};
typedef std::unordered_map<unsigned int, rt_entry> read_tmap;

class ReadProcessor {
public:
  ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, fastaIndex& findex,LKFilter& like_filter,std::stringstream& logfile_track,int thread_num,EarlyTaxonomy& earlytax,bool collapse);
  ~ReadProcessor();
  char *buffer;
  int thread_num;
  size_t bufsize;
  bool paired;
  const MinCollector& tc;
  const KmerIndex& index;
  const fastaIndex& findex;
  const double& illumina_version;
  LKFilter& like_filter;
  EarlyTaxonomy& earlytax;
  MasterProcessor& mp;
  std::stringstream& logfile_track;
  //int numreads;
  std::string out_base;
  double min_logl;
  bool collapse;
  bool fail;
  bool likeplot;
 
  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  //std::vector<std::vector<int>> newEcs;

  //std::vector<int> counts;

  void WriteInt(std::vector<std::unordered_map<double,std::vector<unsigned int> > >& clikes,gzFile& full_outfile,gzFile& single_outfile,std::stringstream& logfile_track,LKFilter& like_filter,std::unordered_map<unsigned int, bool>& c_z_support,gzFile& likefail,std::unordered_map<int, std::string>& read_order,gzFile& maxlikes);
  double Alignment(const char* &s1,const char* &q1,int& ient,const uint64_t& linelength, bool& paired,rt_entry& map_entry,std::string& seqstring);
  double Alignment(const char* &s1,const char* &s2,const char* &q1, const char* &q2,int& ient, const uint64_t& linelength, bool& paired,rt_entry& map_entry,std::string& seqstring);
  std::string getSeqString(int& ient); 
  
  std::unordered_map<double, std::vector<unsigned int > > interpretTmap(read_tmap& tmap);
  void operator()();
  thread_output processBuffer();
  void clear();
};

//double GetMValue(std::vector<double>& likelihoods);

char twinletter(const char* oletter);

#endif // KARP_PROCESSREADS_H
