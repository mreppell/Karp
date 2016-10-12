#ifndef TAXONOMY_HPP
#define TAXONOMY_HPP

#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "karp_like_multi.hpp"
#include "fastaIndex.h"

struct m_node {

  std::string real_name;
  int tax_lvl;
  mutable std::unordered_map<std::string, bool> children;
  mutable double freq;

};

struct collapse_node {
  std::string label;
  std::string id;
};

struct EarlyTaxonomy {
  std::unordered_map<int, unsigned int> label_map;
  std::unordered_map<std::string, unsigned int> name_map;
  std::unordered_map<unsigned int, unsigned int> size_map;
  std::string t_delim;
  int t_format;
  std::unordered_map<std::string, m_node> taxonomy_graph;
  std::unordered_map<int,collapse_node> collapse_map;

  void BuildEarlyTaxonomy(std::vector<std::string>& taxonomy_files,fastaIndex& findex);
  void buildTaxonomyGraph(HLK& likelihoods);
  void output(HLK& likelihoods,std::ofstream& outfile);
  void output2(HLK& likelihoods,std::ofstream& outfile);
  void outputnode(std::string& nodename,std::string rankID,int order,std::ofstream& outfile,double& nummapped);
  void GetDelimeter(std::string& tline);
  std::vector<std::string> TaxName(std::vector<std::string>& taxa);
};

// struct Taxonomy {
  
//   std::unordered_map<std::string, std::vector<std::string> > tax_map;
//   std::string& output_format;
//   std::vector<std::string>& taxonomy_files;
//   
//   HLK& likelihoods;
//   fastaIndex& findex;

//   Taxonomy(ProgramOptions& opt,HLK& likelihoods,fastaIndex& findex) : likelihoods(likelihoods), findex(findex), output_format(opt.format), taxonomy_files(opt.taxonomy_files) {
//     createTaxMap(taxonomy_files,likelihoods,findex);
//   }

//   void createTaxMap(std::vector<std::string>& taxonomy_files,HLK& likelihoods,fastaIndex& findex);
 

// };

typedef std::pair<std::string, m_node> m_entry;
typedef std::pair< std::string, std::vector<std::string> > tax_entry;
typedef std::unordered_map<std::string, m_node >::const_iterator mgraph_get;
typedef std::unordered_map<std::string, bool>::const_iterator mnode_get;
#endif // END TAXONOMY_HPP 
