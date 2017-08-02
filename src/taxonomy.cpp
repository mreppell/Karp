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

void EarlyTaxonomy::readinfo(gzFile& rinfo_file,gzFile& rinfo_out,fastaIndex& findex,std::vector<std::string>& transfasta,gzFile& failfile_in) {
 

  std::stringstream rinfo_ss;
  //rinfo_ss.precision(3);
  rinfo_ss << "#Read Information\n#MapFail = Reads that failed to pseudoalign\n#LikeFail = Reads that failed to pass likelihood filter\n#EM_Unmap = Reads that uniquely pseudoaligned to reference absent from pool after EM convergence\n#Even_PostEM = Reads that offer no evidence for particular references after EM algorithm, and are distributed according to posterior probabilities\n#Using reference IDs from:\n";
  for (int ii=0;ii<transfasta.size();++ii) {
    rinfo_ss << "\t" << transfasta[ii] << "\n";
  }
  rinfo_ss << "#Format Read_Name Probability:ReferenceID\n"; 
  uint64_t buffsize = 1ULL<<20;
  int still_reading = 0;
  char *buff = new char[buffsize];
  int num_read = 0;
  char *amp;
  char *com;
  char *spc;
  char *incurr;
  std::string current = "";

  std::unordered_map<int,std::string> rv_names;
  for (auto it=findex.ref_names.begin();it!=findex.ref_names.end();++it) {
    std::string v1 = it->first;
    int v2 = it->second;
    std::pair<int,std::string> entry(v2,v1);
    rv_names.insert(entry);
  }

  while ((still_reading = gzread(rinfo_file,buff,buffsize))!=0) {
    
    std::string active = std::string(buff);
    std::string new_active = active.substr(0,still_reading);
    active = new_active;
    new_active.clear();
    
    size_t amp = active.find('&');
    size_t curr = 0;
    while (amp!=std::string::npos) {
      if (current=="") {
    	std::string full = active.substr(curr,(amp-curr));
    	curr = amp+1;
    	size_t incurr = full.find(' ');
    	std::string rname = full.substr(0,incurr);
  	rinfo_ss << rname;
  	std::vector<double> percs;
  	std::vector<std::string> na;
  	size_t spc = full.find(' ',incurr+1);
    	while (spc!=std::string::npos) {
  	  std::string part = full.substr(incurr+1,(spc-(incurr+1)));
    	  size_t com = part.find(',');	  
    	  std::string perc = part.substr(0,com);
  	  if (perc!="PF" && perc!="EM") {
  	    std::string ptax = part.substr(com+1,part.size()-(com+1));
  	    int t_ptax = atoi(ptax.c_str());
  	    auto get = rv_names.find(t_ptax);
  	    if (get!=rv_names.end()) {
  	      percs.push_back(atof(perc.c_str()));
  	      na.push_back(get->second);
  	    } else {
  	      std::cerr << "Error finding " << t_ptax << " in fastaIndex for readinfo\n";
  	    }
  	  } else {
  	    std::cerr << "Finding PF/EM for multimapped read: " << rname << std::endl;
  	  }
    	  incurr = spc;
  	  spc = full.find(' ',incurr+1);
    	}
  	std::string opart = full.substr(incurr+1,full.size()-(incurr+1));
  	size_t com = opart.find(',');	  
  	std::string perc = opart.substr(0,com);
  	if (perc!="PF" && perc!="EM") {
  	  std::string ptax = opart.substr(com+1,opart.size()-(com+1));
  	  int t_ptax = atoi(ptax.c_str());
  	  auto get = rv_names.find(t_ptax);
  	  if (get!=rv_names.end()) {
  	    percs.push_back(atof(perc.c_str()));
  	    na.push_back(get->second);
  	  } else {
  	    std::cerr << "Error finding " << t_ptax << " in fastaIndex for readinfo\n";
  	  }
  	  if (percs.size()==1) {
  	    rinfo_ss << " " << percs[0] << ":" << na[0] << "\n";
  	  } else {
  	    double max = 0;
	    std::string max_id = "";
  	    for (int ii=0;ii<percs.size();++ii) {
  	      if (percs[ii] > max) {
  		max = percs[ii];
		max_id = na[ii];
  	      }
  	    }
  	    if (max >= 0.9995) {
	      rinfo_ss << " 1:" << max_id << " <5e-4 total prob:";
	      for (int ii=0;ii<na.size();++ii) {
		if (na[ii]!=max_id) {
		  rinfo_ss << " " << na[ii];
		}
	      }
	      rinfo_ss << "\n";
	    } else {
	      std::sort (percs.begin(),percs.end());
	      std::reverse (percs.begin(),percs.end());
	      for (int ii=0;ii<percs.size();++ii) {
		rinfo_ss << " " << percs[ii] << ":" << na[ii];
	      }
	      rinfo_ss << "\n";
	    }
  	  }
  	} else {
  	  if (perc=="PF") {
  	    rinfo_ss << " EM_Unmap\n";
  	  }
  	  if (perc=="EM") {
  	    rinfo_ss << " Even_PostEM\n";
  	  }
  	}
  	amp = active.find('&',curr+1);
      } else {
    	std::string prefull = active.substr(curr,(amp-curr));
    	std::string full = current + prefull;
    	current = "";
    	curr = amp+1;
    	size_t incurr = full.find(' ');
    	std::string rname = full.substr(0,incurr);
  	rinfo_ss << rname;
  	std::vector<double> percs;
  	std::vector<std::string> na;
  	size_t spc = full.find(' ',incurr+1);
    	while (spc!=std::string::npos) {
  	  std::string part = full.substr(incurr+1,(spc-(incurr+1)));
    	  size_t com = part.find(',');	  
    	  std::string perc = part.substr(0,com);
  	  if (perc!="PF" && perc!="EM") {
  	    std::string ptax = part.substr(com+1,part.size()-(com+1));
  	    int t_ptax = atoi(ptax.c_str());
  	    auto get = rv_names.find(t_ptax);
  	    if (get!=rv_names.end()) {
  	      percs.push_back(atof(perc.c_str()));
  	      na.push_back(get->second);
  	    } else {
  	      std::cerr << "Error finding " << t_ptax << " in fastaIndex for readinfo\n";
  	    }
  	  } else {
  	    std::cerr << "Finding PF/EM for multimapped read: " << rname << std::endl;
  	  }
    	  incurr = spc;
  	  spc = full.find(' ',incurr+1);
    	}
  	std::string opart = full.substr(incurr+1,full.size()-(incurr+1));
  	size_t com = opart.find(',');	  
  	std::string perc = opart.substr(0,com);
  	if (perc!="PF" && perc!="EM") {
  	  std::string ptax = opart.substr(com+1,opart.size()-(com+1));
  	  int t_ptax = atoi(ptax.c_str());
  	  auto get = rv_names.find(t_ptax);
  	  if (get!=rv_names.end()) {
  	    percs.push_back(atof(perc.c_str()));
  	    na.push_back(get->second);
  	  } else {
  	    std::cerr << "Error finding " << t_ptax << " in fastaIndex for readinfo\n";
  	  }
  	  if (percs.size()==1) {
  	    rinfo_ss << " " << percs[0] << ":" << na[0] << "\n";
  	  } else {
  	    double max = 0;
	    std::string max_id = "";
  	    for (int ii=0;ii<percs.size();++ii) {
  	      if (percs[ii] > max) {
  		max = percs[ii];
		max_id = na[ii];
  	      }
  	    }
  	    if (max >= 0.9995) {
	      rinfo_ss << " 1:" << max_id << " <5e-4 total prob:";
	      for (int ii=0;ii<na.size();++ii) {
		if (na[ii]!=max_id) {
		  rinfo_ss << " " << na[ii];
		}
	      }
	      rinfo_ss << "\n";
	    } else {
	      std::sort (percs.begin(),percs.end());
	      std::reverse (percs.begin(),percs.end());
	      for (int ii=0;ii<percs.size();++ii) {
		rinfo_ss << " " << percs[ii] << ":" << na[ii];
	      }
	      rinfo_ss << "\n";
	    }
  	  }
  	} else {
  	  if (perc=="PF") {
  	    rinfo_ss << " EM_Unmap\n";
  	  }
  	  if (perc=="EM") {
  	    rinfo_ss << " Even_PostEM\n";
  	  }
  	}
  	amp = active.find('&',curr+1);
      }
    }
    if (current=="") {
      current = active.substr(curr,active.size()-curr);
    } else {
      std::string precurrent = active.substr(curr,active.size()-curr);
      current = current + precurrent;
    }
  }
  if (current!="") {    
    size_t incurr = current.find(' ');
    std::string rname = current.substr(0,incurr);
    rinfo_ss << rname;
    std::vector<double> percs;
    std::vector<std::string> na;
    size_t spc = current.find(' ',incurr+1);
    while (spc!=std::string::npos) {
      std::string part = current.substr(incurr+1,(spc-(incurr+1)));
      size_t com = part.find(',');	  
      std::string perc = part.substr(0,com);
      if (perc!="PF" && perc!="EM") {
	std::string ptax = part.substr(com+1,part.size()-(com+1));
	int t_ptax = atoi(ptax.c_str());
	auto get = rv_names.find(t_ptax);
	if (get!=rv_names.end()) {
	  percs.push_back(atof(perc.c_str()));
	  na.push_back(get->second);
	} else {
	  std::cerr << "Error finding " << t_ptax << " in fastaIndex for readinfo\n";
	}
      } else {
	std::cerr << "Finding PF/EM for multimapped read: " << rname << std::endl;
      }
      incurr = spc;
      spc = current.find(' ',incurr+1);
    }
  
    std::string opart = current.substr(incurr,current.size()-incurr);
    size_t com = opart.find(',');	  
    std::string perc = opart.substr(0,com);
    std::string ptax = opart.substr(com+1,opart.size()-(com+1));
    if (perc!="PF" && perc!="EM") {
      std::string ptax = opart.substr(com+1,opart.size()-(com+1));
      int t_ptax = atoi(ptax.c_str());
      auto get = rv_names.find(t_ptax);
      if (get!=rv_names.end()) {
	percs.push_back(atof(perc.c_str()));
	na.push_back(get->second);
      } else {
	std::cerr << "Error finding " << t_ptax << " in fastaIndex for readinfo\n";
      }
      if (percs.size()==1) {
	rinfo_ss << " " << percs[0] << ":" << na[0] << "\n";
      } else {
	double max = 0;
	std::string max_id = "";
	for (int ii=0;ii<percs.size();++ii) {
	  if (percs[ii] > max) {
	    max = percs[ii];
	    max_id = na[ii];
	  }
	}
	if (max >= 0.9995) {
	  rinfo_ss << " 1:" << max_id << " <5e-4 total prob:";
	  for (int ii=0;ii<na.size();++ii) {
	    if (na[ii]!=max_id) {
	      rinfo_ss << " " << na[ii];
	    }
	  }
	  rinfo_ss << "\n";
	} else {
	  std::sort (percs.begin(),percs.end());
	  std::reverse(percs.begin(),percs.end());
	  for (int ii=0;ii<percs.size();++ii) {
	    rinfo_ss << " " << percs[ii] << ":" << na[ii];
	  }
	  rinfo_ss << "\n";
	}
      }
    } else {
      if (perc=="PF") {
	rinfo_ss << " EM_Unmap\n";
      }
      if (perc=="EM") {
	rinfo_ss << " Even_PostEM\n";
      }
    }
  }
  delete[] buff;

  still_reading = 0;
  char *abuff = new char[buffsize];
  char *acol;
  char *aspc;
  char *acurr;
  char *nwline;
  char *aincurr;


  while ((still_reading = gzread(failfile_in,abuff,buffsize))!=0) {
    std::string active = std::string(buff);
    std::string new_active = active.substr(0,still_reading);
    active = new_active;
    new_active.clear();
    
    size_t nwline = active.find('\n');
    size_t curr = 0;
    while (nwline!=std::string::npos) {
      std::string prefull = active.substr(curr,nwline-curr);
      std::string full;
      if (current=="") {
	full = prefull;
      } else {
	full = current + prefull;
	current = "";
      }
      size_t acol = full.find(':');
      std::string lid = full.substr(0,acol);
      if (acol!=std::string::npos) {
	std::string reads = full.substr(acol+1,full.size()-(acol+1));
	size_t aincurr = 0;
	size_t aspc = reads.find(' ');
	while (aspc!=std::string::npos) {
	  std::string read = reads.substr(aincurr,aspc-aincurr);
	  rinfo_ss << read;
	  if (lid=="Pseudomapping_Fails") {
	    rinfo_ss << " MapFail\n";
	  }
	  if (lid=="Likelihood_Filter_Fails") {
	    rinfo_ss << " LikeFail\n";
	  }
	  aincurr=aspc+1;
	  aspc = reads.find(' ',aincurr);
	}
	std::string last_read = reads.substr(aincurr,reads.size()-aincurr);
	if (last_read!="") {
	  rinfo_ss << last_read;
	  if (lid=="Pseudomapping_Fails") {
	    rinfo_ss << " MapFail\n";
	  }
	  if (lid=="Likelihood_Filter_Fails") {
	    rinfo_ss << " LikeFail\n";
	  }
	}
      } else {
	std::cerr << "Error, didn't find : in line of fail file\n";
      }
      curr = nwline + 1;
      nwline = active.find('\n',curr);
    }
    if (current=="") {
      current = active.substr(curr,active.size()-curr);
    } else {
      std::string precurrent = active.substr(curr,active.size()-curr);
      current = current + precurrent;
    }
  }
  if (current!="") {
    std::cerr << "Problem with failfile to readinfo, current not empty\n";
  }

  delete[] abuff;
  gzwrite(rinfo_out,rinfo_ss.str().c_str(),rinfo_ss.str().size());
  

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
    if (prec >= 3) {
      prec+=2;
    } else {
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
      if (prec >= 3) {
	prec+=2;
      } else {
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
