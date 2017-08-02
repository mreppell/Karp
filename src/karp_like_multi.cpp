#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <time.h>
#include <limits>
#include <sstream>
#include <cmath>

#include "karp_like_multi.hpp"

std::vector<int> interpretCigar(const char* cigar_string,double& deletions) {

  int num_to_push = 0;
  std::vector<int> matches; 

  for (int ii=0;ii<strlen(cigar_string);++ii) {
    
    if (isdigit(cigar_string[ii])) {
      if (num_to_push==0) {
	num_to_push = cigar_string[ii] - '0';
      } else {
	num_to_push = num_to_push*10;
	num_to_push+=cigar_string[ii] - '0';
      }
    } else {
      if (cigar_string[ii]=='=') {
	for (int jj=0;jj<num_to_push;++jj) {
	  matches.push_back(0);
	}
	num_to_push = 0;
      }
      if (cigar_string[ii]=='S') {
	for (int jj=0;jj<num_to_push;++jj) {
	  matches.push_back(1);
	}
	num_to_push = 0;
      }
      if (cigar_string[ii]=='D') {
	deletions+=num_to_push;
	num_to_push = 0;
      }
      if (cigar_string[ii]=='I' || cigar_string[ii]=='X') {
	for (int jj=0;jj<num_to_push;++jj) {
	  matches.push_back(1);
	}
	num_to_push = 0;
      }
      if (num_to_push!=0) {
	std::cerr << "Error, unrecognized character in Cigar output from aligner: " << cigar_string[ii] << std::endl << "Ignoring bases\n";
	for (int jj=0;jj<num_to_push;++jj) {
	  matches.push_back(2);
	}
	num_to_push = 0;
      }
    }

  }
   
 return matches;
}  

double getError(const char& qual,double mod,std::vector<int>& check) {
  double q = qual - mod;
  if (q < 0) {
    check[0] = 1;
  }
  if (q > 42) {
    check[1] = 1;
  }
  if (mod==64) {
    if (q==0) {
      check[2] = 1;
    }
    if (q==1) {
      check[3] = 1;
    }
    //Legacy from harp, without this q==2 should return ~0.63, but with illumina 1.5+ 'B' represented "Unknown quality score" so given greater error probability than default value
    if (q==2) {
      return 0.75;
    }
  }
  return pow(10, -q/10);
}


double calcLikelihood(std::string& cigar_string, const char* quals,const double& illumina_version) 
 {

   double deletions = 0;
   //Use cigar string from alignment to get matches/mismatches
   std::vector<int> cigar_matches = interpretCigar(cigar_string.c_str(),deletions);
   
   if (cigar_matches.size() != strlen(quals)) {
     std::cerr << "Cigar String does not match base quality length";
     exit(1);
   }

   double mod = illumina_version;
       
   double result = 0;
   const size_t length = cigar_matches.size();
   std::vector<int> check; for (int i=0;i<4;++i) {check.push_back(0);}

   for (size_t i=0; i<length; ++i)
    {
      double p_error = getError(quals[i],mod,check);
      double logp;
      if (cigar_matches[i]==0) {logp = log(1-p_error);}
      if (cigar_matches[i]==1) {logp = log(p_error/3.0);}
      if (cigar_matches[i]==2) {logp = 0;}
      result += logp;
    }

   if (check[3]==1) {
     std::cerr << "Error, detected quality value inconsistent with illumina versions that use Phred+64: A\nCheck version of quality coding being used\n";
   }
   if (check[2]==1) {
     std::cerr << "Error, detected quality value inconsistent with illumina versions that use Phred+64: @\nCheck version of quality coding being used\n";
   }
   if (check[1]==1) {
     std::cerr << "Warning, detected base quality > 42. Verify quality score coding is correct\n";
   }
   if (check[0]==1) {
     std::cerr << "Error, quality scores <0 detected\nPlease enter correct Phred scale factor\n";
     exit(1);
   }

   //Add score for deletions
   //Gap score -2 = same as adding mismatch at site with QC score 9 ~ 40.6% chance of error
   result+= -2*deletions;
   return result;
}


FPTOut HLK::FPTfn(std::unordered_map<unsigned int,double>& current_freqs,bool get_to_min) {

  LikelihoodProcessor MainLikeProcessor(current_freqs,threads,convergence_threshold,num_single,prior_counts,out_base,zero_support,min_logl);
  MainLikeProcessor.ProcessLikes(current_freqs,current_frequency_cutoff);

   if (minimum_frequency_cutoff > 0) {
     if (current_frequency_cutoff + frequency_cutoff_iteration < minimum_frequency_cutoff) {
       if (get_to_min==false) {
	 current_frequency_cutoff+= frequency_cutoff_iteration;
       } else {
	 current_frequency_cutoff = minimum_frequency_cutoff;
       }
     } else {
       if (current_frequency_cutoff!=minimum_frequency_cutoff) {
	 current_frequency_cutoff = minimum_frequency_cutoff;
       }
     }
   }

   FPTOut em_return;
   em_return.freqs = MainLikeProcessor.freqs;
   em_return.final_logl = MainLikeProcessor.final_logl;
   
   if (minimum_frequency_cutoff > 0) {
     if (current_frequency_cutoff==minimum_frequency_cutoff) {
       em_return.at_min = true;
     } else {
       em_return.at_min = false;
     }
   } else {
     em_return.at_min = true;
   }
   return em_return;
   
}

void HLK::singleReadInfo(int& thread,gzFile& catchfile,std::unordered_map<unsigned int, double>& current_freqs) {

  std::stringstream catchinfo;
  std::stringstream curr_sfile_ss;
  curr_sfile_ss << out_base << "_thread" << thread << "_readinfo_singles.gz";
  gzFile curr_sfile = gzopen(curr_sfile_ss.str().c_str(),"rb");
  uint64_t altbuffsize = 1ULL<<20;
  int still_reading = 0;
  char *sbuffer = new char[altbuffsize];
  int num_read = 0;
  char *amp;
  char *pipe; char *curr;;
  std::string current = "";

 while ((still_reading = gzread(curr_sfile,sbuffer,altbuffsize)) > 0) {
    
    std::string active = std::string(sbuffer);
    std::string new_active = active.substr(0,still_reading);
    active = new_active;
    new_active.clear();
    
    size_t pipe = active.find('&');
    size_t curr = 0;
    while (pipe!=std::string::npos) { 
      if (current=="") {
	std::string full = active.substr(curr,(pipe-curr));
	curr = pipe+1;
	size_t amp = full.find('|');
	if (amp!=std::string::npos) {
	  std::string read = full.substr(0,amp);
	  std::string taxa = full.substr(amp+1,full.size());
	  unsigned int  tfind = std::strtoul(taxa.c_str(),NULL,0);
	  auto get = current_freqs.find(tfind);
	  if (get!=current_freqs.end()) {	  
	    catchinfo << read << " 1," << taxa << "&";
	  } else {
	    catchinfo << read << " PF," << taxa << "&";
	  }
	} else {
	  std::cerr << "Error with single file parsing, thread " << thread << ": " << full << std::endl;
	}
	pipe = active.find('&',curr);
	++num_read;
      } else {
	std::string pre_full = active.substr(curr,(pipe-curr));
	std::string full = current + pre_full;
	curr = pipe+1;
	size_t amp = full.find('|');
	if (amp!=std::string::npos) {
	  std::string read = full.substr(0,amp);
	  std::string taxa = full.substr(amp+1,full.size());
	  unsigned int tfind = std::strtoul(taxa.c_str(),NULL,0);
	  auto get = current_freqs.find(tfind);
	  if (get!=current_freqs.end()) {	  
	    catchinfo << read << " 1," << taxa << "&";
	  } else {
	    catchinfo << read << " PF," << taxa << "&";
	  }
	} else {
	  std::cerr << "Error with single file parsing, thread " << thread << ": " << full << std::endl;
	}
	pipe = active.find('&',curr);
	++num_read;
	current = "";
      }
    }
    current = active.substr(curr,active.size()-curr);
  }
  if (current!="") {
    size_t amp = current.find('|',0);    
    if (amp!=std::string::npos) {
      std::string read = current.substr(0,amp);
      std::string taxa = current.substr(amp+1,current.size());
      unsigned int tfind = std::strtoul(taxa.c_str(),NULL,0);
      auto get = current_freqs.find(tfind);
      if (get!=current_freqs.end()) {	  
	catchinfo << read << " 1," << taxa << "&";
      } else {
	catchinfo << read << " PF," << taxa << "&";
      }
    } else {
      std::cerr << "Error with single file parsing, thread " << thread << ": " << current << std::endl;
    }
  }

  gzwrite(catchfile,catchinfo.str().c_str(),catchinfo.str().size());  
  
  gzclose(curr_sfile);
  delete[] sbuffer;
  
}

void HLK::intermedReadInfo(std::unordered_map<unsigned int, double>& current_freqs,std::stringstream& outfile,int& cthread,gzFile& catchfile) {
  
  std::stringstream catchinfo;
  std::stringstream read_info;
  read_info << out_base << "_thread" << cthread << "_readinfo.gz";
  gzFile readfile = gzopen(read_info.str().c_str(),"rb");
  gzFile infile = gzopen(outfile.str().c_str(),"rb");

  std::vector<std::string> thread_reads;
  uint64_t altbuffsize = 1ULL<<20;
  int a_reading = 0;
  char *firstbuffer = new char[altbuffsize];
  int anum_read = 0;
  char *spc; char *curr;
  std::string acurrent = "";
  
  while ((a_reading = gzread(readfile,firstbuffer,altbuffsize)) > 0) {
    
    std::string active = std::string(firstbuffer);
    std::string new_active = active.substr(0,a_reading);
    active = new_active;
    new_active.clear();
    
    size_t spc = active.find(' ');
    size_t curr = 0;
    while (spc!=std::string::npos) { 
      if (acurrent=="") {
	std::string cread = active.substr(curr,(spc-curr));
	curr = spc+1;
	thread_reads.push_back(cread);
	spc = active.find(' ',curr);
	++anum_read;
      } else {
	std::string pre_cread = active.substr(curr,(spc-curr));
	std::string cread = acurrent + pre_cread;
	curr = spc+1;
	thread_reads.push_back(cread);
	acurrent = "";
	spc = active.find(' ',curr);
	++anum_read;
      }
    }
    acurrent = active.substr(curr,active.size()-curr);
  }
  if (acurrent!="") {
    thread_reads.push_back(acurrent);
  }

  gzclose(readfile);
  delete[] firstbuffer;

  char *altbuffer = new char[altbuffsize];

  int still_reading = 0;
  int num_read = 0;
  char *amp; char *pipe;char *at;char *inpipe;char *inat;char *inc1;
  std::string current = "";
    
  while ((still_reading = gzread(infile,altbuffer,altbuffsize)) > 0) {

    std::string active = std::string(altbuffer);
    std::string new_active = active.substr(0,still_reading);
    active = new_active;
    new_active.clear();
    
    size_t amp = active.find('&');
    size_t pipe = active.find('|');
    size_t at = active.find('@');
    size_t c1 = 0;
    std::unordered_map<unsigned int, double> current_read;
    double current_like = 0;
    unsigned int current_id;

    while (amp!=std::string::npos) {

      if (current.compare("")!=0) {
	size_t inpipe = current.find('|');
	size_t inat = current.find('@');
	size_t inc1 = 0;
	  
	while (inpipe!=std::string::npos || inat!=std::string::npos) {
	  if (inpipe < inat && inpipe!=std::string::npos) {
	    current_like = atof(current.substr(inc1,inpipe-inc1).c_str());
	    inc1 = inpipe + 1;
	    inpipe = current.find('|',inc1);
	  }
	  if (inat < inpipe && inat!=std::string::npos) {
	    current_id = atoi(current.substr(inc1,inat-inc1).c_str());
	    std::pair<unsigned int,double> entry(current_id,current_like);
  	     
	    current_read.insert(entry);
	    inc1 = inat + 1;
	    inat = current.find('@',inc1);
	  }
	}
	std::string new_current = current.substr(inc1,current.size()-inc1);
	//std::cout << "New Current: " << new_current << std::endl;
	current = new_current;
      }
	  
      while (pipe < amp || at < amp) {
	  
	if (pipe < at && pipe!=std::string::npos) {
	  if (current.compare("")==0) {
	    current_like = atof(active.substr(c1,pipe-c1).c_str());
	    c1 = pipe + 1;
	    pipe = active.find('|',c1);
	  } else {
	    std::string chelp = current + active.substr(c1,pipe-c1);
	    //std::cout << "chelp: " << chelp << std::endl;
	    current_like = atof(chelp.c_str());
	    c1 = pipe + 1;
	    pipe = active.find('|',c1);
	    current = "";
	  }
	}
	if (at < pipe && at!=std::string::npos) {
	  if (current.compare("")==0) {
	    current_id = atoi(active.substr(c1,at-c1).c_str());
	    std::pair<unsigned int, double> entry(current_id,current_like);
	    //std::cout << "Entry1: " << current_id << " " << current_like << std::endl;
	    current_read.insert(entry);
	    c1 = at + 1;
	    at = active.find('@',c1);
	  } else {
	    std::string chelp = current + active.substr(c1,at-c1);
	    current_id = atoi(chelp.c_str());
	    std::pair<unsigned int, double> entry(current_id,current_like);
	    //std::cout << "Entry2: " << current_id << " " << current_like << std::endl;
	    current_read.insert(entry);
	    c1 = at + 1;
	    at = active.find('@',c1);
	    current = "";
	  }
	}
      }
       
      double marginal = base_marginal;
      
     for (auto it=current_read.begin();it!=current_read.end();++it) {
	auto get = current_freqs.find(it->first);
	if (get!=current_freqs.end()) {
	  double val = exp(it->second)*(get->second);
	  marginal = marginal + val - exp(min_logl)*(get->second);
	}
      }

      catchinfo << thread_reads[num_read];
      bool see_something = false;
      for (auto it=current_read.begin();it!=current_read.end();++it) {
	auto get = current_freqs.find(it->first);
	if (get!=current_freqs.end()) {
	  see_something = true;
	  double val = exp(it->second)*(get->second);
	  double rval = val/marginal;
	  std::cout.precision(3);
	  catchinfo << " " << rval << "," << get->first;
	}
      }
      if (see_something==false) {
	catchinfo << " EM,None";
      }
      catchinfo << '&'; 

      current_read.clear();
      c1 = amp + 1;
      amp = active.find('&',c1);
      num_read++;
    }
    if (current.compare("")==0) {
      current = active.substr(c1,active.size()-c1);
    } else {
      std::string new_current = current + active.substr(0,active.size());
      current = new_current;
    }

  }
  gzwrite(catchfile,catchinfo.str().c_str(),catchinfo.str().size());

  gzclose(infile);
  delete[] altbuffer;
    
}

void LikelihoodProcessor::ProcessLikes(std::unordered_map<unsigned int,double>& current_freqs,double& frequency_cutoff) {
  
  double base_marginal = 0;
  for (auto get=current_freqs.begin();get!=current_freqs.end();++get) {
    //std::cout << get->first << "(" << get->second << ") ";
    base_marginal+= exp(min_logl)*(get->second);
  }
  //std::cout << std::endl;
  //std::cout << "Base Marginal: " << base_marginal << std::endl;

  std::vector<std::thread> workers;
  for (int i = 0; i < threads; i++) {  
    workers.emplace_back(std::thread(GetAbundance,this,i,current_freqs,base_marginal,out_base,min_logl));    
  }
  
  // let the workers do their thing
  for (int i = 0; i < threads; i++) {
    workers[i].join(); //wait for them to finish
  }

  double total_reads = 0;
  for (int i=0;i < threads;++i) {
    total_reads+=read_counts[i];
  }
  
  for (auto it=current_freqs.begin();it!=current_freqs.end();++it) {
    auto get = prior_counts.find(it->first);
    if (get!=prior_counts.end()) {
      total_reads+=get->second;
    }
  }

  std::unordered_map<unsigned int,double> new_freqs;
  //double sumsum2 = 0;  

  for (auto it=current_freqs.begin();it!=current_freqs.end();++it) {
    double ab_count = 0;    
    //std::cout << it->first << "|| ";
    for (int i=0;i<threads;++i) {
      auto get = abundance_counts[i].find(it->first);
      if (get==abundance_counts[i].end()) {
	std::cerr << "Error in getting counts from each thread, shouldn't happen\n";
	exit(1);
      }
      //std::cout << "(" << i << ")" << get->second << " "; 
      ab_count+=get->second;
    }
    auto got = prior_counts.find(it->first);
    if (got!=prior_counts.end()) {
      ab_count+=got->second;
    }
    //std::cout << ab_count/total_reads << std::endl;
    std::pair<unsigned int,double> entry(it->first,(ab_count/total_reads));
    //sumsum2+=ab_count/total_reads;
    // if (it->first==171) {
    //   std::cout << "Ab_count " << ab_count << " total_reads " << total_reads << " ";
    //   double val2 = ab_count/total_reads;
    //   std::cout << it->first << " " << val2 << std::endl;
    // }
    new_freqs.insert(entry);
  }
  //std::cout << "New Freq Sum: " << sumsum2 << " Total Reads " << total_reads << std::endl;
 
  if (frequency_cutoff > 0) {
    cutoff_low_frequencies(current_freqs,new_freqs,frequency_cutoff);
  }
   
  double final_like = 0;
  for (int i=0;i<threads;++i) {
    final_like+=marginal_likes[i];
  }
  
  freqs = new_freqs;
  final_logl = final_like;

}


SquaremOutput HLK::SquareEM(std::unordered_map<unsigned int, double>& current_freqs){
  double res,parnorm,kres;
  FPTOut pcpp,p1cpp,p2cpp,pnew,ptmp;
  std::vector<double> q1,q2,sr2,sq2,sv2,srv;
  double sr2_scalar,sq2_scalar,sv2_scalar,srv_scalar,alpha,stepmin,stepmax;

  int iter,feval;
  bool conv,extrap,get_to_min;
  stepmin=EMcontrol.stepmin0;
  stepmax=EMcontrol.stepmax0;
  
  iter=1;
  pcpp.freqs = current_freqs;
  pnew.freqs = current_freqs;
  feval=0;
  conv=true;
  get_to_min = false;

  const long int parvectorlength0=pcpp.freqs.size();

  while(feval<EMcontrol.maxiter){
    //Step 1
    extrap = true;
    try{p1cpp=FPTfn(pcpp.freqs,get_to_min);feval++;}
    catch(...){
      std::cerr<<"Error in EM function evaluation";
      exit(1);
    }

    sr2_scalar=0;
    for (auto get1=pcpp.freqs.begin();get1!=pcpp.freqs.end();++get1) {
      auto get2 = p1cpp.freqs.find(get1->first);
      if (get2!=p1cpp.freqs.end()) {
  	sr2_scalar+=pow((get2->second)-(get1->second),2);
      } else {
  	sr2_scalar+=pow(get1->second,2);
      }
    }

    //std::cout << "PCPP " << pcpp.freqs.size() << " P1CPP " << p1cpp.freqs.size);
    
    if (sqrt(sr2_scalar)<EMcontrol.tol) {
      	pcpp = p1cpp;
      if (p1cpp.at_min==true) {
	break;
      } else {
	get_to_min = true;
      }
    }

    //Step 2
    try{p2cpp=FPTfn(p1cpp.freqs,get_to_min);feval++;}
    catch(...){
      std::cout<<"Error in EM function evaluation";
      exit(1);
    }
    sq2_scalar=0;
    for (auto get1=p1cpp.freqs.begin();get1!=p1cpp.freqs.end();++get1) {
      auto get2 = p2cpp.freqs.find(get1->first);
      if (get2!=p2cpp.freqs.end()) {
  	sq2_scalar+=pow((get2->second)-(get1->second),2);
      } else {
  	sq2_scalar+=pow(get1->second,2);
      }
    }
    //std::cout << " P2CPP " << p2cpp.freqs.size();
    if (sq2_scalar<EMcontrol.tol) {
      pcpp = p2cpp;
      if (p2cpp.at_min==true) {
	break;
      } else {
	get_to_min = true;
      }
    }
     
    sv2_scalar=0;
    srv_scalar=0;
    for (auto get1=pcpp.freqs.begin();get1!=pcpp.freqs.end();++get1) {
      auto get2 = p1cpp.freqs.find(get1->first);
      auto get3 = p2cpp.freqs.find(get1->first);
      double g2 = 0;double g3 = 0;
      if (get2!=p1cpp.freqs.end()) {g2 = get2->second;}
      if (get3!=p2cpp.freqs.end()) {g3 = get3->second;}
      sv2_scalar+=pow(g3-2*g2+(get1->second),2);
      srv_scalar+=(g3-2*g2+(get1->second))*(g2-(get1->second));
    }
     
    //Step 3 Proposing new value
    switch (EMcontrol.method){
    case 1: alpha= -srv_scalar/sv2_scalar;
    case 2: alpha= -sr2_scalar/srv_scalar;
    case 3: alpha= sqrt(sr2_scalar/sv2_scalar);
    }
        
    alpha=std::max(stepmin,std::min(stepmax,alpha)); 
    double new_weight = 0;
    //std::cout << "CFC " << current_frequency_cutoff << " Beginning update size " << pnew.freqs.size();
    for (auto get1=pcpp.freqs.begin();get1!=pcpp.freqs.end();++get1) {
      auto get2 = p1cpp.freqs.find(get1->first);
      auto get3 = p2cpp.freqs.find(get1->first);
      auto newget = pnew.freqs.find(get1->first);
      double g2 = 0;double g3 = 0; 
      if (get2!=p1cpp.freqs.end()) {g2 = get2->second;}
      if (get3!=p2cpp.freqs.end()) {g3 = get3->second;}
      double newfreq = (get1->second)+2.0*alpha*(g2-(get1->second))+pow(alpha,2)*(g3-2*g2+(get1->second));
      if (newfreq > current_frequency_cutoff) {
	if (newget!=pnew.freqs.end()) {
	  newget->second = newfreq;
	} else {
	  double newfreqentry = newfreq;
	  std::pair<unsigned int, double> nfe(get1->first,newfreqentry);
	  pnew.freqs.insert(nfe);
	}
	new_weight+=newfreq;
      } else {
	if (newget!=pnew.freqs.end()) {
	  pnew.freqs.erase(newget);
	}
      }
    }
    //std::cout << " PNEW " << pnew.freqs.size() << std::endl;
    

    //std::cout << " Ending: " << pnew.freqs.size() << " weight " << new_weight;
    //double freq_sum = 0;
    for (auto it1=pnew.freqs.begin();it1!=pnew.freqs.end();++it1) {
      it1->second = it1->second/new_weight;
      //freq_sum+=it1->second;
    }
    // std::cout << " New FreqSum: " << freq_sum << std::endl;
    
    //Step 4 stabilization
    if(std::abs(alpha-1)>0.01){
      try{ptmp=FPTfn(pnew.freqs,get_to_min);feval++;}
      catch(...){
  	pnew=p2cpp;
  	if(alpha==stepmax){
  	  stepmax=std::max(EMcontrol.stepmax0,stepmax/EMcontrol.mstep);
  	}
  	alpha=1;
  	extrap=false;
  	if(alpha==stepmax){stepmax=EMcontrol.mstep*stepmax;}
  	if(stepmin<0 && alpha==stepmin){stepmin=EMcontrol.mstep*stepmin;}
  	pcpp=pnew;
  	if(EMcontrol.trace){std::cerr<<"Residual: "<<res<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
  	iter++;
  	continue;//next round in while loop
      }
      res=0;
      for (auto get1=pnew.freqs.begin();get1!=pnew.freqs.end();++get1) {
  	auto get2 = ptmp.freqs.find(get1->first);
  	if (get2!=ptmp.freqs.end()) {
  	  res+=pow((get2->second)-(get1->second),2);
  	} else {
  	  res+=pow(get1->second,2);
  	}
      }

      res=sqrt(res);
      parnorm=0;
      for (auto get1=p2cpp.freqs.begin();get1!=p2cpp.freqs.end();++get1) {
  	parnorm+=pow((get1->second),2);
      }
      parnorm=sqrt(parnorm/p2cpp.freqs.size());
      kres=EMcontrol.kr*(1+parnorm)+sq2_scalar;
      if(res <= kres){
  	pnew=ptmp;
      }else{
  	pnew=p2cpp;
  	if(alpha==stepmax){stepmax=EMcontrol.mstep*stepmax;}
  	alpha=1;
  	extrap=false;
      }
    }
       	    
    if(alpha==stepmax){stepmax=EMcontrol.mstep*stepmax;}
    if(stepmin<0 && alpha==stepmin){stepmin=EMcontrol.mstep*stepmin;}

    pcpp=pnew;
    if(EMcontrol.trace){std::cout<<"Residual: "<<res<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
    iter++;
  }

  if (feval >= EMcontrol.maxiter){conv=false;}

  SquaremOutput output;
  output.freqs = pcpp.freqs;
  output.iter = iter;
  output.pfevals = feval;
  output.convergence = conv;
  output.final_logl = pcpp.final_logl;

  return output;
}

void HLK::estimate_haplotype_frequencies() {

  std::unordered_map<unsigned int,double> current_freqs;

  GetPriors(logfile_track);

  double ifreq = initializeFreqs(current_freqs,zero_support,prior_counts,nummapped);

  //  for (auto it = current_freqs.begin();it!=current_freqs.end();++it) {
  //   std::cout << "(" << it->first << "):" << it->second << " ";
  // }
  // std::cout << std::endl;
    
  frequency_cutoff_iteration = minimum_frequency_cutoff/50.0;
  current_frequency_cutoff = frequency_cutoff_iteration;
  
  SquaremOutput em_return = SquareEM(current_freqs);

  if (em_return.convergence==true) {
    std::cerr << "EM converged in " << em_return.pfevals << " with final log(like) = " << em_return.final_logl << std::endl;
    logfile_track << "EM converged in " << em_return.pfevals << " with final log(like) = " << em_return.final_logl << std::endl;
    freqs = em_return.freqs;
  } else {
    std::cerr << "EM failed to converge in maximum number of iterations(" << em_return.pfevals << "), using final iteration frequencies, but these may be inaccurate" << std::endl;
    logfile_track << "EM failed to converge in maximum number of iterations(" << em_return.pfevals << "), using final iteration frequencies, but these may be inaccurate" << std::endl;
    freqs = em_return.freqs;
  }
  
  for (auto it=prior_counts.begin();it!=prior_counts.end();++it) {
    auto get = freqs.find(it->first);
    if (get==freqs.end()) {
      num_single_fail+=it->second;
    }
  }
  nummapped-=num_single_fail;

  std::stringstream readinfo_catch_all;
  readinfo_catch_all << out_base << "_readinfo_intermediate.txt.gz";
  gzFile catchfile = gzopen(readinfo_catch_all.str().c_str(),"wb");

  for (int ii=0;ii<threads;++ii) {
    std::stringstream outfilename;
    std::stringstream singleoutfile;
    outfilename << out_base << "_thread" << ii << ".lks";
    singleoutfile << out_base << "_thread" << ii << "_singles.lks";
    if (readinfo) {
      intermedReadInfo(freqs,outfilename,ii,catchfile);
      singleReadInfo(ii,catchfile,freqs);
    }
    std::stringstream ri_file1;
    std::stringstream ri_file2;
    ri_file1 << out_base << "_thread" << ii << "_readinfo.gz";
    ri_file2 << out_base << "_thread" << ii << "_readinfo_singles.gz";
    if ( remove(ri_file1.str().c_str()) != 0) {
      perror("Error deleting readinfo file");
    }
    if ( remove(ri_file2.str().c_str()) != 0) {
      perror("Error deleting readinfo file");
    }
    if ( remove(outfilename.str().c_str()) != 0) {
      perror( "Error deleting file");
    } 
    if ( remove(singleoutfile.str().c_str()) != 0) {
      perror( "Error deleting file");
    } 

  }
  gzclose(catchfile);
  if (!readinfo) {
    if ( remove(readinfo_catch_all.str().c_str()) != 0) {
      perror ("Error deleting empty readinfo intermediate file");
    }
  }
   
  std::cerr << "After EM " << num_single_fail << " reads uniquely mapped to a reference taxa with a frequency below the minimum threshold and were not assigned\n" << nummapped << " reads remain and were assigned a taxonomic label" << std::endl;
  logfile_track << "After EM " << num_single_fail << " reads uniquely mapped to a reference taxa with a frequency below the minimum threshold and were not assigned\n" << nummapped << " reads remain and were assigned a taxonomic label" << std::endl;
}

double HLK::initializeFreqs(std::unordered_map<unsigned int,double>& current_freqs,std::unordered_map<unsigned int,bool>& zero_support,std::unordered_map<unsigned int,double>& prior_counts,double& nummapped) {
   
  double denom = 0;
  for (auto i=zero_support.begin();i!=zero_support.end();++i) {
    if (i->second==true) {
      //     std::cout << "Post: " << i->first << " " << i->second << std::endl;
      denom++;
    }
  }
  //std::cout << "Denom: " << denom << std::endl;

  for (auto i=zero_support.begin();i!=zero_support.end();++i) {
    //std::cout << i->first << std::endl;
    if (i->second==true) {
      auto inget = prior_counts.find(i->first);
      if (inget!=prior_counts.end()) {
	double pfreq = (inget->second/nummapped) + remaining_freq/denom;
	std::pair<unsigned int, double> entr(i->first,pfreq);
	current_freqs.insert(entr);
      } else {
	double pfreq = remaining_freq/denom;
	std::pair<unsigned int, double> entr(i->first,pfreq);
	current_freqs.insert(entr);
      }
    }
  }
  return remaining_freq/denom;
}


void  cutoff_low_frequencies(std::unordered_map<unsigned int,double>& f_old,std::unordered_map<unsigned int,double>& f_new,double& min_freq) {

  double f_new_sum = 0;
  //std::cout << "Before there were: " << f_new.size() << " refs"<< std::endl;
  for (auto it=f_old.begin();it!=f_old.end();++it) {
    auto get = f_new.find(it->first);
    if (get==f_new.end()) {
      std::cerr << "Error, in frequency cutoff step failed to find reference " << it->first << std::endl;
    }
    if (get->second < min_freq) {
      f_new.erase(get);
    } else {
      f_new_sum+= get->second;
    }
  }
  // std::cerr << std::endl << std::endl;
  for (auto it=f_new.begin();it!=f_new.end();++it) {
    it->second = (it->second)/f_new_sum;
  }
  
  //std::cerr << "From " << f_old.size() << " -> "<< f_new.size() << " Cutoff " << min_freq << std::endl;
}
/////////Idea is to create modified version of function below, call on final results frequencies with read names


void GetAbundance(LikelihoodProcessor* LP,int thread_num,std::unordered_map<unsigned int,double> current_freqs,double base_marginal,std::string out_base,double min_logl) {

  double log_sum_marginal = 0;
  std::stringstream outfilename;
  outfilename << out_base << "_thread" << thread_num << ".lks";
  gzFile infile = gzopen(outfilename.str().c_str(),"rb");
  if (!infile) {
    std::cerr << "Unable to open non-single likelihood file for thread " << thread_num << std::endl;
    exit(1);
  }

  uint64_t altbuffsize = 1ULL<<20;
  char *altbuffer = new char[altbuffsize];

  int still_reading = 0;
  int num_read = 0;
  char *amp; char *pipe;char *at;char *inpipe;char *inat;char *inc1;
  std::string current = "";
    
  while ((still_reading = gzread(infile,altbuffer,altbuffsize)) > 0) {

    std::string active = std::string(altbuffer);
    std::string new_active = active.substr(0,still_reading);
    active = new_active;
    new_active.clear();
    
     size_t amp = active.find('&');
     size_t pipe = active.find('|');
     size_t at = active.find('@');
     size_t c1 = 0;
     std::unordered_map<unsigned int, double> current_read;
     double current_like = 0;
     unsigned int current_id;

     while (amp!=std::string::npos) {

       if (current.compare("")!=0) {
	 size_t inpipe = current.find('|');
	 size_t inat = current.find('@');
	 size_t inc1 = 0;
	  
	 while (inpipe!=std::string::npos || inat!=std::string::npos) {
	   if (inpipe < inat && inpipe!=std::string::npos) {
	     current_like = atof(current.substr(inc1,inpipe-inc1).c_str());
	     inc1 = inpipe + 1;
	     inpipe = current.find('|',inc1);
	   }
	   if (inat < inpipe && inat!=std::string::npos) {
	     current_id = atoi(current.substr(inc1,inat-inc1).c_str());
	     std::pair<unsigned int,double> entry(current_id,current_like);
	     //std::cout << "Repeat from current: " << current_id << " " << current_like << std::endl;
	     current_read.insert(entry);
	     inc1 = inat + 1;
	     inat = current.find('@',inc1);
	   }
	 }
	 std::string new_current = current.substr(inc1,current.size()-inc1);
	 //std::cout << "New Current: " << new_current << std::endl;
	 current = new_current;
       }
	  
       while (pipe < amp || at < amp) {
	  
	 if (pipe < at && pipe!=std::string::npos) {
	   if (current.compare("")==0) {
	     current_like = atof(active.substr(c1,pipe-c1).c_str());
	     c1 = pipe + 1;
	     pipe = active.find('|',c1);
	   } else {
	     std::string chelp = current + active.substr(c1,pipe-c1);
	     //std::cout << "chelp: " << chelp << std::endl;
	     current_like = atof(chelp.c_str());
	     c1 = pipe + 1;
	     pipe = active.find('|',c1);
	     current = "";
	   }
	 }
	 if (at < pipe && at!=std::string::npos) {
	   if (current.compare("")==0) {
	     current_id = atoi(active.substr(c1,at-c1).c_str());
	     std::pair<unsigned int, double> entry(current_id,current_like);
	     //std::cout << "Entry1: " << current_id << " " << current_like << std::endl;
	     current_read.insert(entry);
	     c1 = at + 1;
	     at = active.find('@',c1);
	   } else {
	     std::string chelp = current + active.substr(c1,at-c1);
	     current_id = atoi(chelp.c_str());
	     std::pair<unsigned int, double> entry(current_id,current_like);
	     //std::cout << "Entry2: " << current_id << " " << current_like << std::endl;
	     current_read.insert(entry);
	     c1 = at + 1;
	     at = active.find('@',c1);
	     current = "";
	   }
	 }
       }
       
       double marginal = base_marginal;
       for (auto it=current_read.begin();it!=current_read.end();++it) {
	 auto get = current_freqs.find(it->first);
	 if (get!=current_freqs.end()) {
	   double val = exp(it->second)*(get->second);
	   marginal = marginal + val - exp(min_logl)*(get->second);

	   //if (it->first==946) {
	   //  std::cout << "Th: " << thread_num << " val " << val << " exp(" << it->second << ") " << get->second << " " << marginal << std::endl;
	   // }
	 }
       }
       log_sum_marginal+=log(marginal);
       
       //if (num_read==284) {
	 //std::cout << "Read " << num_read << " ";
         //	for ( auto get = current_read.begin();get!=current_read.end();++get) {
         //	  std::cout << get->first << " " << get->second << std::endl;
         //	}
         //}

       double sumsum = 0;
       for (auto it=LP->abundance_counts[thread_num].begin();it!=LP->abundance_counts[thread_num].end();++it) {
       	 auto get=current_freqs.find(it->first);
       	 if (get==current_freqs.end()) {
       	   std::cerr << "Didn't find abundance in current freqs, this shouldn't be happening\n";
       	 }
       	 auto got=current_read.find(it->first);
       	 if (got==current_read.end()) {
       	   double val = (exp(min_logl)*(get->second))/marginal;
       	   it->second+= val;
       	   sumsum+=val;
       	 } else {
       	   double val = (exp(got->second)*(get->second))/marginal;
       	   it->second+=val;
	   sumsum+=val;
       	 }
       }
       //std::cout << std::endl;
       //std::cout << "Thread " << thread_num << " read " << num_read << " sum: " << sumsum << " marginal: " << marginal << std::endl; 

       current_read.clear();
       c1 = amp + 1;
       amp = active.find('&',c1);
       num_read++;
     }
     if (current.compare("")==0) {
       current = active.substr(c1,active.size()-c1);
     } else {
       std::string new_current = current + active.substr(0,active.size());
       current = new_current;
     }
     //std::cout << "Current at end: " << current << std::endl;

  }
  LP->read_counts[thread_num] = num_read;
  //std::cout << "LSM Thread " << thread_num << ": " << log_sum_marginal << std::endl;
  LP->marginal_likes[thread_num] = log_sum_marginal;

  gzclose(infile);
  delete[] altbuffer;

}
	


// std::unordered_map<unsigned int,double> HLK::GetPriors(std::stringstream& logfile_track) {
void HLK::GetPriors(std::stringstream& logfile_track) {
   
  //std::unordered_map<unsigned int,double> prior_counts;
 
   uint64_t buff1size = 1ULL<<20;
   //uint64_t buffsize = 128;
   char *buffer1 = new char[buff1size];
 
  for (int th=0;th<threads;++th) {
    
    std::stringstream single_outfilename;
    single_outfilename << out_base << "_thread" << th << "_singles.lks";
    gzFile infile = gzopen(single_outfilename.str().c_str(),"rb");
    if (!infile) {
      std::cerr << "Unable to open intermediate singles file for thread " << th << std::endl;
      exit(1);
    }
     
    int still_reading = 0;
    int num_read = 0;
    char *amp; char *pipe;char *at;char *inpipe;char *inat;char *inc1;
    std::string current = "";

    while ((still_reading = gzread(infile,buffer1,buff1size)) > 0) {
      
      std::string active = std::string(buffer1);      
      std::string new_active = active.substr(0,still_reading);
      active = new_active;
      new_active.clear();

      size_t amp = active.find('&');
      size_t c1 = 0;
      unsigned int current_id;

      while (amp!=std::string::npos) {

	 if (current.compare("")==0) {
	   current_id = atoi(active.substr(c1,amp-c1).c_str());
	   //std::cout << current_id << std::endl;
	   auto get = prior_counts.find(current_id);
	   if (get==prior_counts.end()) {
	     std::pair<unsigned int,double> entry(current_id,1.0);
	     prior_counts.insert(entry);
	   } else {
	     get->second++;
	   }
	 } else {
	   std::string chelp = current + active.substr(c1,amp-c1).c_str();
	   //std::cout << "chelp " << chelp << std::endl;
	   current_id = atoi(chelp.c_str());
	   //std::cout << current_id << std::endl;
	   current = "";
	   auto get = prior_counts.find(current_id);
	   if (get==prior_counts.end()) {
	     std::pair<unsigned int,double> entry(current_id,1.0);
	     prior_counts.insert(entry);
	   } else {
	     get->second++;
	   }
	 }
	 ++num_single;
	 c1 = amp+1;
	 amp = active.find('&',c1);
      }

      if (current.compare("")==0) {
	current = active.substr(c1,active.size()-c1);
      } else {
	std::string new_current = current + active.substr(0,active.size());
	current = new_current;
      }

      //std::cout << "Current " << current << std::endl;

    }
  }

  std::cerr << num_single << " reads uniquely mapped to a reference and were assigned before EM algorithm\n";
  logfile_track << num_single << " reads uniquely mapped to a reference and were assigned before EM algorithm\n";
  delete[] buffer1;

  // double total_pfreq = 0;
  // for (auto get=prior_counts.begin();get!=prior_counts.end();++get) {
  //   double cfreq = get->second/nummapped;
  //   total_pfreq+=cfreq;
  //   std::pair<unsigned int,double> ee(get->first,cfreq);
  //   prior_freqs.insert(ee);
  // }
  remaining_freq = 1.0 - num_single/nummapped;
  
}
