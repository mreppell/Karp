#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "taxonomy.hpp"
#include "common_k.h"

void errorMessage() {
  std::cerr << "Error, unsupported taxonomy file format\n";
}

std::string  gettlev(int& x) {

  if (x==7) {return "species";}
  if (x==6) {return "genus";}
  if (x==5) {return "family";}
  if (x==4) {return "order";}
  if (x==3) {return "class";}
  if (x==2) {return "phylum";}
  if (x==1) {return "kingdom";}

  std::cerr << "Unkown taxnomic level " << x << std::endl;
  exit(1);
}

std::vector<std::string> EarlyTaxonomy::TaxName(std::vector<std::string>& taxa) {
      
  std::string in_delim = "__";
  std::vector<std::string> new_taxa;
  int still_good = 0;
  for (int ii=0;ii<taxa.size();++ii) {
    if (!taxa[ii].empty()) {
      if (t_format==0) {
	std::vector<std::string> f0;
	split(taxa[ii],in_delim,f0);
	if (f0.size()>1) {
	  if (still_good==0) {	        
	    new_taxa.push_back(f0[1]);
	  } else {
	    std::cerr << "Error " << taxa[ii] << " seen after higher taxa empty\n";
	  }
	} else {
	  ++still_good;
	}
      } else {
	new_taxa.push_back(taxa[ii]);
      }
    } else {
      break;
    }
  }
  taxa = new_taxa;
  if(taxa.size()>=7) {
    return taxa;
  } else {
    int lowest_t_level = taxa.size() + 1;
    std::string lowest_t_unit = taxa[taxa.size()-1];
    while (lowest_t_level<7) {
      ++lowest_t_level;
      std::string tlev = gettlev(lowest_t_level);
      std::string new_label = "unknown_" + lowest_t_unit + "_" + tlev;
      taxa.push_back(new_label);
    }	
    return taxa;
  }  
}

void EarlyTaxonomy::GetDelimeter(std::string& tline) {
  int format = 0;
  std::string delim2 = "; ";
  std::string delim3 = ";";
  std::string delim4 = "__";
  t_delim = delim2;
  
  std::vector<std::string> first_line_taxa;
  split(tline,delim2,first_line_taxa);
   
  if (first_line_taxa.size()==1) {
    std::vector<std::string> fl2;
    split(tline,delim3,fl2);
    if (fl2.size()==1) {errorMessage();exit(1);}
   
    //split(tline,delim3,first_line_taxa2);
    //if (first_line_taxa2.size()==1) {errorMessage();exit(1);}
    std::vector<std::string> u_check;
    split(fl2[0],delim4,u_check);
    if (u_check.size()>2) {errorMessage();exit(1);}
    if (u_check.size()<2) {format=1;}
    t_delim = delim3;
  } else {
    std::vector<std::string> u_check;
    split(first_line_taxa[0],delim4,u_check);
    if (u_check.size()>2) {errorMessage();exit(1);}
    if (u_check.size()<2) {format=1;}
  }      
  t_format = format;
}

void EarlyTaxonomy::buildTaxonomyGraph(HLK& likelihoods) {

  m_node root_node;
  std::string nroot = "Root";
  //std::string uk = "Root__Unknown";
  root_node.real_name = nroot;
  root_node.tax_lvl = 0;
  root_node.freq = 1.0;
  m_entry root(nroot,root_node);

  taxonomy_graph.insert(root);

  int unique_addition = 1;
  
  for (auto it=name_map.begin(); it!=name_map.end();++it) {
    unsigned int id_num = it->second;
    std::string name = it->first;
    //std::cout << it->first << std::endl;

    auto getter = likelihoods.freqs.find(id_num);
    if (getter!=likelihoods.freqs.end()) {
     
      std::vector<std::string> taxa;
      split(name,t_delim,taxa);
      std::vector<std::string> tax;
      if (taxa.size()>1) {
	tax = TaxName(taxa);
      } else {
	tax.push_back(taxa[0]);
      }
      int c_tax = tax.size();
      std::stringstream uniquename;
      for (int cc=0;cc<c_tax;++cc) {
	if (cc==0) {
	  //std::stringstream uniquename;
	  uniquename << tax[cc];
	  mgraph_get get1 = taxonomy_graph.find("Root");
	  mnode_get get2 = (get1->second).children.find(uniquename.str());
	    
	  if (get2 == (get1->second).children.end()) {
	    m_node c_entry;
	    c_entry.real_name = tax[cc];
	    c_entry.tax_lvl = 1;
	    c_entry.freq = getter->second;

	    m_entry node(uniquename.str(),c_entry);
	    taxonomy_graph.insert(node);
	    std::pair<std::string, bool> new_child(uniquename.str(),true);
	    (get1->second).children.insert(new_child);
	  } else {
	    mgraph_get get3 = taxonomy_graph.find(uniquename.str());
	    (get3->second).freq+=getter->second;
	  }
	} else {
	  //std::stringstream unique_parentname = uniquename;
	  //unique_parentname << tax[cc-1] << cc;
	  mgraph_get get1 = taxonomy_graph.find(uniquename.str());
	  if (get1 == taxonomy_graph.end()) {std::cerr << "Something wrong with taxonomy graph construction, this shouldn't be happening!!!\n";exit(1);}
	  //std::stringstream uniquename;
	  uniquename << "__" << tax[cc];
	  mnode_get get2 = (get1->second).children.find(uniquename.str());
	  if (get2 == (get1->second).children.end()) {
	    m_node c_entry;
	    c_entry.real_name = tax[cc];
	    c_entry.tax_lvl = cc + 1;
	    c_entry.freq = getter->second;
	    m_entry node(uniquename.str(),c_entry);
	    taxonomy_graph.insert(node);
	    std::pair<std::string, bool> new_child(uniquename.str(),true);
	    (get1->second).children.insert(new_child);
	  } else {
	    mgraph_get get3 = taxonomy_graph.find(uniquename.str());
	    (get3->second).freq+=getter->second;
	  }
	}
      }
    }
  }
}
  
void EarlyTaxonomy::outputnode(std::string& nodename,std::string rankID,int order,std::ofstream& outfile,double& nummapped) {

  
  //std::cout << nodename << std::endl;
  mgraph_get get1 = taxonomy_graph.find(nodename);
 
  std::vector<std::string> unknown_check;
  std::string delim = "_";
  std::string check_name = get1->second.real_name;
  split(check_name,delim,unknown_check);
  //std::cerr << check_name << " " << unknown_check[0]; // << std::endl;
  if (unknown_check[0].compare("unknown")!=0) {
    

    std::stringstream rID;
    double ecount = (get1->second).freq*nummapped;
    int prec = (int) log10(ecount);
    if (prec > 4) {
      prec+=1;
    } else {
      if (prec==3) {
	prec+=2;
      }
      if (prec==2) {
	prec+=3;
      }
      if (prec < 2) {
	prec = 4;
      }
    }

    if (order==0) {
      outfile << (get1->second).tax_lvl << "\t" << order << "\t" << (get1->second).real_name << "\t" << (get1->second).children.size() << "\t" << ecount << std::endl;    
      rID << order;
    } else {
      rID << rankID << "." << order;
      
      outfile << (get1->second).tax_lvl << "\t" << rID.str() << "\t" << (get1->second).real_name << "\t" << (get1->second).children.size() << "\t" << std::setprecision(prec) << ecount << std::endl;
      
    }
    
    int kid_count = 1;
    for (auto get2 = (get1->second).children.begin();get2 != (get1->second).children.end();++get2) {
      std::string inname = get2->first;
      std::string inRD = rID.str();
      outputnode(inname,inRD,kid_count,outfile,nummapped);
      ++kid_count;
    }
  }
  
}

void EarlyTaxonomy::output2(HLK& likelihoods,std::ofstream& outfile) {
  
  outfile << "Label\tExpectedCounts\tTaxa\n";
  for (auto get=likelihoods.freqs.begin();get!=likelihoods.freqs.end();++get) {
    int old_fient = (int) get->first;
    double ecount = (get->second)*likelihoods.nummapped;
    auto got = collapse_map.find(old_fient);
    if (got==collapse_map.end()) {
      std::cerr << "Error finding " << old_fient << " in taxonomy, skipping\n";
    } else {
      int prec = (int) log10(ecount);
      if (prec > 4) {
	prec+=1;
      } else {
	if (prec==3) {
	  prec+=2;
	}
	if (prec==2) {
	  prec+=3;
	}
	if (prec < 2) {
	  prec = 4;
	}
      }  
      outfile << (got->second).id << "\t" << std::setprecision(prec) << ecount << "\t" << (got->second).label << std::endl;
    }
  }
}


void EarlyTaxonomy::output(HLK& likelihoods,std::ofstream& outfile) {
   
    buildTaxonomyGraph(likelihoods);

    outfile << "taxlevel\trankID\ttaxon\tdaughterlevels\texpected_counts\n";
    std::string r1 = "Root";
    std::string rID1 = "0";
    outputnode(r1,rID1,0,outfile,likelihoods.nummapped);
    
} 
  
void EarlyTaxonomy::BuildEarlyTaxonomy(std::vector<std::string>& taxonomy_files,fastaIndex& findex) {

  std::string delim1 = "\t";
  for (int file=0;file<taxonomy_files.size();++file) {
    const char* fn = taxonomy_files[file].c_str();
    std::ifstream myfile (fn);
    std::string line;
    
    unsigned int cat_num = 1;
    bool got_delimeter = false;

    if (myfile.is_open()) {

      while(getline(myfile,line)) {
	std::vector<std::string> first_line;
	split(line,delim1,first_line);
	if (first_line.size() != 2) {errorMessage();exit(1);}
	if (got_delimeter==false) {
	  GetDelimeter(first_line[1]);
	  got_delimeter=true;
	}

	std::unordered_map<std::string, int >::const_iterator got1 = findex.ref_names.find(first_line[0]);
	if (got1!=findex.ref_names.end()) {	
	  int f_ient = got1->second;
	  std::unordered_map<std::string, unsigned int>::const_iterator get2 = name_map.find(first_line[1]);
	  if (get2!=name_map.end()) {
	    unsigned int category = get2->second;
	    auto oget = size_map.find(category);
	    oget->second++;
	    std::pair<int, unsigned int> entry(f_ient,category);
	    label_map.insert(entry);

	    collapse_node coll_entry;
	    coll_entry.label = get2->first;
	    coll_entry.id = got1->first;
	    std::pair<int, collapse_node > collpair(f_ient,coll_entry);
	    collapse_map.insert(collpair);


	  } else {
	  
	    std::pair<std::string, unsigned int> newlab(first_line[1],cat_num);
	    name_map.insert(newlab);
	    std::pair<unsigned int, unsigned int> nentry(cat_num,1);
	    size_map.insert(nentry);
	    std::pair<int, unsigned int> entry(f_ient,cat_num);
	    label_map.insert(entry);
	    ++cat_num;

	    collapse_node coll_entry;
	    coll_entry.label = first_line[1];
	    coll_entry.id = got1->first;
	    std::pair<int, collapse_node > collpair(f_ient,coll_entry);
	    collapse_map.insert(collpair);

	  }
	}
      }
      myfile.close();
    } else {
      std::cerr << "Unable to open taxonomy file " << taxonomy_files[file] << std::endl;
      exit(1);
    }
  }
}
