#ifndef KARP_MINCOLLECTOR_H
#define KARP_MINCOLLECTOR_H

#include "common_k.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

#include "KmerIndex.h"
//#include "weights.h"

const int MAX_FRAG_LEN = 1000;

struct MinCollector {

MinCollector(KmerIndex& ind, const ProgramOptions& opt)
:
  index(ind),
    counts(index.ecmap.size(), 0),
    min_range(opt.min_range), resc_thresh(opt.resc_thresh), lenient(opt.multi_match),
    k(opt.k)
  {}
  
  int collect(std::vector<std::pair<KmerEntry,int>>& v1,
              std::vector<std::pair<KmerEntry,int>>& v2,
              bool nonpaired=false);

  int collect(std::vector<std::pair<KmerEntry,int>>& v1) {
    std::vector<std::pair<KmerEntry,int>> dummy;
    return collect(v1,dummy,true);

  }
  int increaseCount(const std::vector<int>& u);
  int decreaseCount(const int ec);


  std::vector<int> printECs(std::vector<std::pair<KmerEntry,int>>& v) const;
  std::vector<int> intersectECs(std::vector<std::pair<KmerEntry,int>>& v) const;
  std::vector<int> altIntersectECs(std::vector<std::pair<KmerEntry,int>>& v) const;
  int intersectKmers(std::vector<std::pair<KmerEntry,int>>& v1,
                    std::vector<std::pair<KmerEntry,int>>& v2, bool nonpaired, std::vector<int> &u) const;
  int altIntersectKmers(std::vector<std::pair<KmerEntry,int>>& v1,
			std::vector<std::pair<KmerEntry,int>>& v2, bool nonpaired, std::vector<int> &u) const;
  
  int findEC(const std::vector<int>& u) const;

  void write(std::ostream& o) {
    for (int id = 0; id < counts.size(); id++) {
      o << id << "\t" << counts[id] << "\n";
    }
  }
  void loadCounts(ProgramOptions& opt);

  KmerIndex& index;
  std::vector<int> counts;
  int min_range;
  const int k;
  int resc_thresh;
  bool lenient;

};

std::vector<int> intersect(const std::vector<int>& x, const std::vector<int>& y);

int hexamerToInt(const char *s, bool revcomp);

#endif // KARP_MINCOLLECTOR_H
