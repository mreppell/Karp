
#include "tabulate.hpp"

void ResTable::tabulate(std::string& all_files) {

  std::string delim = ",";
  std::string delim2 = "\t";

  const char* bifn = all_files.c_str();
  std::ifstream filelist (bifn);
  if (!filelist.is_open()) {
    std::cerr << "Unable to open file with list of samples to tabulate: " << all_files << std::endl;
    exit(1);
  }
  std::string outline;
  std::vector<std::string> names;
  std::vector<std::string> sfiles;

  while(getline(filelist,outline)) {
    std::vector<std::string> outfiles;
    split(outline,delim2,outfiles);
    if (outfiles.size()!=2) {
      outfiles.empty();
      split(outline,delim,outfiles);
      if (outfiles.size()!=2) {
	std::cerr << "Incorrect format of sample list. Should be two tab or comma separated columns\nObserved: " << outline << std::endl;
	exit(1);
      }
    }
    names.push_back(outfiles[0]);
    sfiles.push_back(outfiles[1]);
  }
  filelist.close();
    
  for (int ii=0;ii<sfiles.size();++ii) {
    
    const char* fn = sfiles[ii].c_str();
    std::ifstream infile (fn);
    
    std::string line;

    double file_totcount = 0;

    if (!infile.is_open()) {
      std::cerr << "Unable to open sample file " << sfiles[ii] << " for tabulation" << std::endl;
      continue;
    } else {
      std::vector<std::string> head_line;
      getline(infile,line);
      split(line,delim2,head_line);
      if (head_line[0].compare("Label")) {
	std::cerr << "Error, incorrect sample file format. Unrecognized column label: "<< head_line[0] << "\n";
	exit(1);
      }
      while(getline(infile,line)) {
	std::vector<std::string> line_vals;
	split(line,delim2,line_vals);
	if (line_vals.size()!=3) {
	  std::cerr << "Error with line: " << line << " unrecognized format\n";
	  continue;
	}
	file_totcount+=atof(line_vals[1].c_str());
      }
      //infile.close();
    }
    infile.clear();
    infile.seekg(0,std::ios::beg);

    std::unordered_map<unsigned int, double> counts;
    double file_mincount = file_totcount*min_f;

    if (!infile.is_open()) {
      std::cerr << "Unable to open sample file " << sfiles[ii] << " for tabulation" << std::endl;
      continue;
    } else {
      std::vector<std::string> head_line;
      getline(infile,line);
      split(line,delim2,head_line);
      while(getline(infile,line)) {
	std::vector<std::string> line_vals;
	split(line,delim2,line_vals);
	if (line_vals.size()!=3) {
	  std::cerr << "Error with line: " << line << " unrecognized format\n";
	  continue;
	}
	if (atof(line_vals[1].c_str()) >= file_mincount) {

	  unsigned int tlab = strtoul(line_vals[0].c_str(), NULL, 0);
	  std::pair<unsigned int,double> entry(tlab,atof(line_vals[1].c_str()));
	  counts.insert(entry);
	  auto get = tax_table.find(tlab);
	  if (get==tax_table.end()) {
	    std::pair<unsigned int,std::string> t_entry(tlab,line_vals[2]);
	    tax_table.insert(t_entry);
	  } else {
	    if (get->second.compare(line_vals[2])!=0) {
	      std::cerr << "Error taxonomy label mismatch ID:" << get->first << " " << get->second << " ne " << line_vals[2] << std::endl;
	    }
	  }
	}
      }
      infile.close();
    }
    std::pair<std::string, std::unordered_map<unsigned int,double> > sample_entry(names[ii],counts);
    otu_table.insert(sample_entry);
		
  }

}

void ResTable::output() {
  
  std::unordered_map<std::string, double> total_counts;
  double total_otus = 0;
  int tot_samples = 0;
  std::vector<double> a_counts;
  double a_count_total = 0;

  //Output table with aggregate raw counts

  std::string countout = out_base + ".combined_counts.txt";
  std::ofstream outcount;
  outcount.open(countout.c_str(), std::ios::out);
  if (!outcount.is_open()) {
    std::cerr << "Unable to open output count file\n";
    exit(1);
  }

  for (auto it=tax_table.begin();it!=tax_table.end();++it) {
    if (it==tax_table.begin()) {
      outcount << it->first;
    } else {
      outcount << "\t" << it->first;
    }
    a_counts.push_back(0);
    total_otus++;
  }
  outcount << std::endl;

  for (auto it=tax_table.begin();it!=tax_table.end();++it) {
    if (it==tax_table.begin()) {
      outcount << it->second;
    } else {
      outcount << "\t" << it->second;
    }
  }
  outcount << std::endl;

  for (auto sam=otu_table.begin();sam!=otu_table.end();++sam) {
    outcount << sam->first;
    double total = 0;
    ++tot_samples;
    int index = 0;
    for (auto it=tax_table.begin();it!=tax_table.end();++it) {
      auto sit = (sam->second).find(it->first);
      if (sit==(sam->second).end()) {
	outcount << "\t0";
      } else {
	outcount << "\t" << sit->second;
	total+=sit->second;
	a_counts[index]+=sit->second;
      }
      ++index;
    }
    outcount << std::endl;
    std::pair<std::string, double> c_entry(sam->first, total);
    total_counts.insert(c_entry);
    a_count_total+=total;
  }
  outcount.close();


  //Calculate summary statistics and output
  std::string stfile = out_base + ".summary_stats.txt";
  std::ofstream statout;
  statout.open(stfile.c_str(), std::ios::out);
  if (!statout.is_open()) {
    std::cerr << "Unable to open output summary statistic file\n";
    exit(1);
  }

  statout << "#########INDV STATISTICS########\n";
  
  std::vector<double> shannon;
  std::vector<double> ginisimpson;
  std::vector<unsigned int> otus;
  for (auto sam=otu_table.begin();sam!=otu_table.end();++sam) {
    double entropy = 0;
    double simpson = 0;
    unsigned int obs = 0;
    auto tot = total_counts.find(sam->first);
    for (auto tax=(sam->second).begin();tax!=(sam->second).end();++tax) {
      double prop = (tax->second)/(tot->second);
      entropy+=prop*log(prop);
      simpson+=prop*prop;
      ++obs;
    }
    otus.push_back(obs);
    shannon.push_back(exp(entropy));
    ginisimpson.push_back(simpson);
  }

 statout << "SAMPLE";
  for (auto it=otu_table.begin();it!=otu_table.end();++it) {
    statout << "\t" << it->first;
  }
  statout << std::endl;
  statout << "D0_DIVERSITY(RICHNESS):";
 for (int ii=0;ii<tot_samples;++ii) {
    statout << "\t" << otus[ii];
  }
  statout <<std::endl;
  statout << "D1_DIVERSITY(SHANNON):";
  for (int ii=0;ii<tot_samples;++ii) {
    statout << "\t" << shannon[ii];
  }
  statout <<std::endl;
  statout << "D2_DIVERSITY(SIMPSON):";
  for (int ii=0;ii<tot_samples;++ii) {
    statout << "\t" << ginisimpson[ii];
  }

  double dgamma1 = 0;
  double dgamma2 = 0;
  for (int ii=0;ii<a_counts.size();++ii) {
    double prop = a_counts[ii]/a_count_total;
    dgamma1 += prop*log(prop);
    dgamma2 += prop*prop;
  }
  dgamma1 = exp(dgamma1);
  
  double d1alpha = 0;
  double d2alpha = 0;
  int alpha_index = 0;
  for (auto it=otu_table.begin();it!=otu_table.end();++it) {
    auto get = total_counts.find(it->first);
    double sprop = get->second/a_count_total;
    d1alpha+=sprop*shannon[alpha_index];
    d2alpha+=sprop*ginisimpson[alpha_index];
    ++alpha_index;
  }
  statout << std::endl;
  
  statout << "\n#########GROUP STATISTICS########\n";
  
  statout << "D1_GAMMA_DIVERSITY: " << dgamma1 << std::endl;
  statout << "D2_GAMMA_DIVERSITY: " << dgamma2 << std::endl;
  statout << "D1_BETA_DIVERSITY: " << dgamma1/d1alpha << std::endl;
  statout << "D2_BETA_DIVERSITY: " << dgamma2/d2alpha << std::endl;

  statout << "\n#########PAIRWISE STATISTICS########\n";
  statout << "PAIRWISE D2_BETA_DIVERSITY:\n";

  std::vector<std::string> name_orders;
  for (auto s1=otu_table.begin();s1!=otu_table.end();++s1) {
    name_orders.push_back(s1->first);
  }
  for (int it1=(name_orders.size()-1);it1>0;--it1) {
    statout << "\t" << name_orders[it1];
  }
  statout << std::endl;

  alpha_index = 0;
  
  for (int ps1=0;ps1<(name_orders.size()-1);++ps1) {
    auto s1 = otu_table.find(name_orders[ps1]);
    if (s1==otu_table.end()) {
      std::cerr << name_orders[ps1] << " Not found! Error!\n";
      exit(1);
    }
    auto t1 = total_counts.find(s1->first);
    int in_index = alpha_index + 1;
    statout << s1->first;

    for (int ps2=(name_orders.size()-1);ps2>ps1;ps2--) {
      auto s2 = otu_table.find(name_orders[ps2]);
      if (s2==otu_table.end()) {
	std::cerr << name_orders[ps2] << " Not found! Error!\n";
	exit(1);
      }
      auto t2 = total_counts.find(s2->first);
      
      std::unordered_map<unsigned int, double> combo_counts;      
      for (auto tx1=(s1->second).begin();tx1!=(s1->second).end();++tx1) {
	std::pair<unsigned int, double> cc_e(tx1->first,tx1->second);
	combo_counts.insert(cc_e);
      }
      for (auto tx2=(s2->second).begin();tx2!=(s2->second).end();++tx2) {
	auto gg = combo_counts.find(tx2->first);
	if (gg==combo_counts.end()) {
	  std::pair<unsigned int, double> cc_e(tx2->first,tx2->second);
	  combo_counts.insert(cc_e);
	} else {
	  gg->second+=tx2->second;
	}
      }

      double d2pairwisegamma = 0;
      for (auto cc1=combo_counts.begin();cc1!=combo_counts.end();++cc1) {
	double prop = (cc1->second)/(t1->second + t2->second);
	d2pairwisegamma += prop*prop;
      }
      
      double d2pairwisebeta = d2pairwisegamma/((t1->second/(t1->second + t2->second))*ginisimpson[ps1] + (t2->second/(t1->second + t2->second))*ginisimpson[ps2]);
      statout << "\t" << d2pairwisebeta;
      ++in_index;
    }
    
    statout << std::endl;
    ++alpha_index;
  }    
  

  statout << "\nPAIRWISE BRAY-CURTIS DISSIMILARITIES:\n";
 
  for (int it1=(name_orders.size()-1);it1>0;--it1) {
    statout << "\t" << name_orders[it1];
  }
  statout << std::endl;
  for (int ps1=0;ps1<(name_orders.size()-1);++ps1) {
    auto s1 = otu_table.find(name_orders[ps1]);
    auto t1 = total_counts.find(s1->first);
    statout << s1->first;
    for (int ps2=(name_orders.size()-1);ps2>ps1;ps2--) {
      auto s2 = otu_table.find(name_orders[ps2]);
      auto t2 = total_counts.find(s2->first);
      double bc_val = 0;
      for (auto tx=tax_table.begin();tx!=tax_table.end();++tx) {
	auto s1tx=(s1->second).find(tx->first);
	auto s2tx=(s2->second).find(tx->first);
	double v1 = 0;
	double v2 = 0;
	if (s1tx!=(s1->second).end()) {
	  v1 = s1tx->second;
	}
	if (s2tx!=(s2->second).end()) {
	  v2 = s2tx->second;
	}
	bc_val+=fabs(v1 - v2);
      }
      bc_val = bc_val / (t1->second + t2->second);
      statout << "\t" << bc_val;
    }
    
    statout << std::endl;
  }
  

  statout.close();

}
