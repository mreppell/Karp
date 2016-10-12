#include "MinCollector.h"
#include <algorithm>

// utility functions

std::vector<int> intersect(const std::vector<int>& x, const std::vector<int>& y) {
  std::vector<int> v;
  auto a = x.begin();
  auto b = y.begin();
  while (a != x.end() && b != y.end()) {
    if (*a < *b) {
      ++a;
    } else if (*b < *a) {
      ++b;
    } else {
      v.push_back(*a);
      ++a;
      ++b;
    }
  }
  return v;
}

int MinCollector::intersectKmers(std::vector<std::pair<KmerEntry,int>>& v1,
                          std::vector<std::pair<KmerEntry,int>>& v2, bool nonpaired, std::vector<int> &u) const {
  std::vector<int> u1 = intersectECs(v1);
  std::vector<int> u2 = intersectECs(v2);

  if (u1.empty() && u2.empty()) {
    return -1;
  }

  // non-strict intersection.
  if (u1.empty()) {
    if (v1.empty()) {
      u = u2;
    } else {
      return -1;
    }
  } else if (u2.empty()) {
    if (v2.empty()) {
      u = u1;
    } else {
      return -1;
    }
  } else {
    u = intersect(u1,u2);
  }

  if (u.empty()) {
    return -1;
  }
  return 1;
}

int MinCollector::collect(std::vector<std::pair<KmerEntry,int>>& v1,
                          std::vector<std::pair<KmerEntry,int>>& v2, bool nonpaired) {
  std::vector<int> u;
  int r = intersectKmers(v1, v2, nonpaired, u);
  if (r != -1) {
    return increaseCount(u);
  } else {
    return -1;
  }
}

int MinCollector::findEC(const std::vector<int>& u) const {
  if (u.empty()) {
    return -1;
  }
  if (u.size() == 1) {
    return u[0];
  }
  auto search = index.ecmapinv.find(u);
  if (search != index.ecmapinv.end()) {
    return search ->second;
  } else {
    return -1;
  }
}

int MinCollector::increaseCount(const std::vector<int>& u) {
  int ec = findEC(u);

  if (u.empty()) {
    return -1;
  } else {
    if (ec != -1) {
      ++counts[ec];
      return ec;
    } else {
      auto necs = counts.size();
      //index.ecmap.insert({necs,u});
      index.ecmap.push_back(u);
      index.ecmapinv.insert({u,necs});
      counts.push_back(1);
      return necs;
    }
  }

}

int MinCollector::decreaseCount(const int ec) {
  assert(ec >= 0 && ec <= index.ecmap.size());
  --counts[ec];
  return ec;
}

struct ComparePairsBySecond {
  bool operator()(std::pair<KmerEntry,int> a, std::pair<KmerEntry,int> b) {
    return a.second < b.second;
  }
};


std::vector<int> MinCollector::intersectECs(std::vector<std::pair<KmerEntry,int>>& v) const {
  if (v.empty()) {
    return {};
  }
  sort(v.begin(), v.end(), [&](std::pair<KmerEntry, int> a, std::pair<KmerEntry, int> b)
       {
         if (a.first.contig==b.first.contig) {
           return a.second < b.second;
         } else {
           return a.first.contig < b.first.contig;
         }
       }); // sort by contig, and then first position


  int ec = index.dbGraph.ecs[v[0].first.contig];
  int lastEC = ec;
  std::vector<int> u = index.ecmap[ec];

  for (int i = 1; i < v.size(); i++) {
    if (v[i].first.contig != v[i-1].first.contig) {
      ec = index.dbGraph.ecs[v[i].first.contig];
      if (ec != lastEC) {
        u = index.intersect(ec, u);
        lastEC = ec;
        if (u.empty()) {
          return u;
        }
      }
    }
  }

  // find the range of support
  int minpos = std::numeric_limits<int>::max();
  int maxpos = 0;

  for (auto& x : v) {
    minpos = std::min(minpos, x.second);
    maxpos = std::max(maxpos, x.second);
  }

  if ((maxpos-minpos + k) < min_range) {
    return {};
  }

  return u;
}

int hexamerToInt(const char *s, bool revcomp) {
  int hex = 0;
  if (!revcomp) {
    for (int i = 0; i < 6; i++) {
      hex <<= 2;
      switch (*(s+i) & 0xDF) {
      case 'A': break;
      case 'C': hex += 1; break;
      case 'G': hex += 2; break;
      case 'T': hex += 3; break;
      default: return -1;
      }
    }
  } else {
    for (int i = 0; i < 6; i++) {
      switch (*(s+i) & 0xDF) {
      case 'A': hex += 3 << (2*i);break;
      case 'C': hex += 2 << (2*i); break;
      case 'G': hex += 1 << (2*i); break;
      case 'T': break;
      default: return -1;
      }
    }
  }
  return hex;
}

int MinCollector::altIntersectKmers(std::vector<std::pair<KmerEntry,int>>& v1,
				    std::vector<std::pair<KmerEntry,int>>& v2, bool nonpaired, std::vector<int> &u) const {

  std::vector<int> u1;
  std::vector<int> u2;

  u1 = altIntersectECs(v1);
  u2 = altIntersectECs(v2);
   
  if (u1.empty()) {
    if (u2.empty()) {
      return -1;
    } else {
      u = u2;
    }
  }
  else if (u2.empty()) {
    u = u1;
  } else {
    u = intersect(u1,u2);
  }

  //If no common references, but few enough that union is below rescue threshold, return union of references
  if (u.empty()) {
    if (u1.size() + u2.size() < resc_thresh && u1.size() + u2.size() > 0) {
      for (int ii=0;ii<u1.size();++ii) {u.push_back(u1[ii]);}
      for (int jj=0;jj<u2.size();++jj) {u.push_back(u2[jj]);}
    } else {
      return -1;
    }
  }
  return 1;
}


std::vector<int> MinCollector::altIntersectECs(std::vector<std::pair<KmerEntry,int>>& v) const {
  if (v.empty()) {
    return {};
  }
  sort(v.begin(), v.end(), [&](std::pair<KmerEntry, int> a, std::pair<KmerEntry, int> b)
       {
         if (a.first.contig==b.first.contig) {
           return a.second < b.second;
         } else {
           return a.first.contig < b.first.contig;
         }
       }); // sort by contig, and then first position
  

  int max_seen = 1;
  std::unordered_map<int, int> seen_multiple_times; 
  std::vector<int> u;

  int ec = index.dbGraph.ecs[v[0].first.contig];
  int change = 0;
  int lastEC = ec;
  std::vector<int> pre_u = index.ecmap[ec];
  //std::cerr << "Read\n(0)";
  for (int ii=0;ii<pre_u.size();++ii) {
    std::pair<int, int> entry(pre_u[ii],1);
    //std::cerr << " " << pre_u[ii];
    seen_multiple_times.insert(entry);
  }
  //std::cerr << std::endl;

  for (int i = 1; i < v.size(); i++) {
    if (v[i].first.contig != v[i-1].first.contig) {
      ec = index.dbGraph.ecs[v[i].first.contig];
      if (ec != lastEC) {
	++change;
	pre_u = index.ecmap[ec];
	//	std::cerr << "(" << i << ")";
	for (int ii=0;ii<pre_u.size();++ii) {
	  //std::cerr << " " << pre_u[ii];
	  auto get = seen_multiple_times.find(pre_u[ii]);
	  if (get!=seen_multiple_times.end()) {
	    (get->second)++;
	    if ( (get->second) > max_seen) {
	      max_seen = get->second;
	    }
	  } else {
	    std::pair<int, int> entry(pre_u[ii],1);
	    seen_multiple_times.insert(entry);
	  }
	  
	 
	}
	//std::cerr << std::endl;
	lastEC = ec;
      }
      
    }
  }
  //std::cerr << "U:";
  
  if (change==0) {
    for (int it=0;it<pre_u.size();++it) {
      //  std::cerr << " " << pre_u[it];
      u.push_back(pre_u[it]);
    
    }
  } else {

    if (max_seen>1) {   
      for (auto it=seen_multiple_times.begin();it!=seen_multiple_times.end();++it) {
	if (it->second==max_seen) {
	  u.push_back(it->first);
	//	std::cerr << " " << it->first;
	}
      }
    }
  }

  // find the range of support
  int minpos = std::numeric_limits<int>::max();
  int maxpos = 0;

  for (auto& x : v) {
    minpos = std::min(minpos, x.second);
    maxpos = std::max(maxpos, x.second);
  }

  if ((maxpos-minpos + k) < min_range) {
    return {};
  }

  return u;
}

// std::vector<int> MinCollector::altIntersectECs(std::vector<std::pair<KmerEntry,int>>& v) const {
//   if (v.empty()) {
//     return {};
//   }
//   sort(v.begin(), v.end(), [&](std::pair<KmerEntry, int> a, std::pair<KmerEntry, int> b)
//        {
//          if (a.first.contig==b.first.contig) {
//            return a.second < b.second;
//          } else {
//            return a.first.contig < b.first.contig;
//          }
//        }); // sort by contig, and then first position


//   std::unordered_map<int, int> seen_multiple_times; 
//   std::vector<int> u;

//   int ec = index.dbGraph.ecs[v[0].first.contig];
//   int change = 0;
//   int lastEC = ec;
//   std::vector<int> pre_u = index.ecmap[ec];
//   //std::cerr << "Read\n(0)";
//   for (int ii=0;ii<pre_u.size();++ii) {
//     std::pair<int, int> entry(pre_u[ii],1);
//     //std::cerr << " " << pre_u[ii];
//     seen_multiple_times.insert(entry);
//   }
//   //std::cerr << std::endl;

//   for (int i = 1; i < v.size(); i++) {
//     if (v[i].first.contig != v[i-1].first.contig) {
//       ec = index.dbGraph.ecs[v[i].first.contig];
//       if (ec != lastEC) {
// 	++change;
// 	pre_u = index.ecmap[ec];
// 	//	std::cerr << "(" << i << ")";
// 	for (int ii=0;ii<pre_u.size();++ii) {
// 	  //std::cerr << " " << pre_u[ii];
// 	  auto get = seen_multiple_times.find(pre_u[ii]);
// 	  if (get!=seen_multiple_times.end()) {
// 	    (get->second)++;
// 	  } else {
// 	    std::pair<int, int> entry(pre_u[ii],1);
// 	    seen_multiple_times.insert(entry);
// 	  }
	  
	 
// 	}
// 	//std::cerr << std::endl;
// 	lastEC = ec;
//       }
      
//     }
//   }
//   //std::cerr << "U:";
  
//   if (change==0) {
//     for (int it=0;it<pre_u.size();++it) {
//       //  std::cerr << " " << pre_u[it];
//       u.push_back(pre_u[it]);
    
//     }
//   } else {

//     for (auto it=seen_multiple_times.begin();it!=seen_multiple_times.end();++it) {
//       if (it->second > 1) {
// 	u.push_back(it->first);
// 	//	std::cerr << " " << it->first;
//       }
//     }
//   }
//   //  std::cerr << std::endl;

//   // find the range of support
//   int minpos = std::numeric_limits<int>::max();
//   int maxpos = 0;

//   for (auto& x : v) {
//     minpos = std::min(minpos, x.second);
//     maxpos = std::max(maxpos, x.second);
//   }

//   if ((maxpos-minpos + k) < min_range) {
//     return {};
//   }

//   return u;
// }

