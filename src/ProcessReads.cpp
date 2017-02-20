/*
#include <zlib.h>
#include "kseq.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>

#include <iostream>
#include <fstream>
#include "MinCollector.h"
#include "common_k.h"
*/
#include <algorithm>
#include <iterator>
#include "fastaIndex.h"
#include "ProcessReads.h"
#include "kseq.h"
#include "ssw_cpp.h"
#include "Kmer.hpp"
#include "karp_like_multi.hpp"
#include <chrono>


char twinletter(const char& oletter) {
  if (oletter=='A') { return 'T';}
  if (oletter=='T') { return 'A';}
  if (oletter=='G') { return 'C';}
  if (oletter=='C') { return 'G';}
  return 'N';
}

void printVector(const std::vector<int>& v, std::ostream& o) {
  o << "[";
  int i = 0;
  for (auto x : v) {
    if (i>0) {
      o << ", ";
    }
    o << x;
    i++;
  }
  o << "]";
}

bool isSubset(const std::vector<int>& x, const std::vector<int>& y) {
  auto a = x.begin();
  auto b = y.begin();
  while (a != x.end() && b != y.end()) {
    if (*a < *b) {
      return false;
    } else if (*b < *a) {
      ++b;
    } else {
      ++a;
      ++b;
    }
  }
  return (a==x.end());
}

uint32_t ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc, fastaIndex& findex,std::stringstream& logfile_track,LKFilter& like_filter,EarlyTaxonomy& earlytax,HLK& likelihoods) {

  // int limit = 131072;
  //char *buf = new char[limit];
  //std::vector<std::pair<const char*, int>> seqs;
  //seqs.reserve(limit/50);

  // need to receive an index map
  std::ios_base::sync_with_stdio(false);

  uint32_t numreads = 0;
  uint32_t numstrict = 0;
  uint32_t numresc = 0;

  bool paired = opt.single_end;

  if (paired) {
    std::cerr << "Running in paired-end mode" << std::endl;
    logfile_track << "Running in paired-end mode" << std::endl;
  } else {
    std::cerr << "Running in single-end mode" << std::endl;
    logfile_track << "Running in single-end mode" << std::endl;
  }

  for (int i = 0; i < opt.files.size(); i += (paired) ? 2 : 1) {
    if (paired) {
      std::cerr   << "Process pair " << (i/2 +1) << "\n\t"  << opt.files[i]
		  << "\n\t" << opt.files[i+1] << std::endl << std::endl; 
      logfile_track  << "Process pair " << (i/2 +1) << "\n\t"  << opt.files[i]
                 << "\n\t" << opt.files[i+1] << std::endl << std::endl; 
    } else {
      std::cerr << "Process file " << i+1 << ": " << opt.files[i] << std::endl;
      logfile_track << "Process file " << i+1 << ": " << opt.files[i] << std::endl;
    }
  }

  // for each file
  std::cerr << "Finding pseudoalignments for the reads ..."; std::cerr.flush();
  logfile_track << "Finding pseudoalignments for the reads ...";
 
  MasterProcessor MP(index, opt, tc, findex, like_filter,logfile_track,earlytax);
  MP.processReads();

  for (int ii=0;ii<opt.threads;++ii) {
    numreads+=MP.numreads[ii];
    numstrict+=MP.numstrict[ii];
    numresc+=MP.numresc[ii];
  }
  uint32_t nummapped = numstrict + numresc;
  uint32_t nleft;

  for (int ii=0;ii<opt.threads;++ii) {
    for (auto it=MP.zsupport[ii].begin();it!=MP.zsupport[ii].end();++it) {
      auto get=likelihoods.zero_support.find(it->first);

      if (get==likelihoods.zero_support.end()) {


	std::pair<unsigned int, bool> nent(it->first,it->second);
	likelihoods.zero_support.insert(nent);
      } else {
	if (get->second==false) {
	  if (it->second==true) {
	    get->second=true;
	  }
	}
      }
    }
  }

  if (opt.fail) {
      std::stringstream pre_full_fail_file;
      pre_full_fail_file << opt.out << ".failed_reads.gz";
      gzFile full_fail_file = gzopen(pre_full_fail_file.str().c_str(),"wb");
      
      std::string row1("Psuedomapping_Fails:");
      std::string row2("\nLikelihood_Filter_Fails:");

      gzwrite(full_fail_file,row1.c_str(),row1.size());
      uint64_t failbuffsize = 1ULL<<20;
      char *failbuff = new char[failbuffsize];
	
      for (int ii=0;ii<opt.threads;++ii) {
	std::stringstream current_fail_file;
	current_fail_file << opt.out << "_thread" << ii << "_mapfails.gz";
	gzFile infile = gzopen(current_fail_file.str().c_str(),"rb");
	if (!infile) {
	  std::cerr << "Unable to open fail file for thread " << ii << std::endl;
	  exit(1);
	}
	int still_reading = 0;
	while ((still_reading = gzread(infile,failbuff,failbuffsize)) > 0 ) {
	  std::string active = std::string(failbuff);
	  std::string new_active = active.substr(0,still_reading);
	  gzwrite(full_fail_file,new_active.c_str(),new_active.size());
	}
	gzclose(infile);
	if ( remove(current_fail_file.str().c_str()) != 0) {
	  perror( "Error deleting mapfail file");
	} 
      }
      
      gzwrite(full_fail_file,row2.c_str(),row2.size());
      for (int ii=0;ii<opt.threads;++ii) {
	std::stringstream current_fail_file;
	current_fail_file << opt.out << "_thread" << ii << "_likefails.gz";
	gzFile infile = gzopen(current_fail_file.str().c_str(),"rb");
	if (!infile) {
	  std::cerr << "Unable to open fail file for thread " << ii << std::endl;
	  exit(1);
	}
	int still_reading = 0;
	while ((still_reading = gzread(infile,failbuff,failbuffsize)) > 0 ) {
	  std::string active = std::string(failbuff);
	  std::string new_active = active.substr(0,still_reading);
	  gzwrite(full_fail_file,new_active.c_str(),new_active.size());
	}
	gzclose(infile);
	if ( remove(current_fail_file.str().c_str()) != 0) {
	  std::cerr << current_fail_file.str().c_str() << " "; 
	  perror("Error deleting likefail file");
	} 
      }
      delete[] failbuff;
      gzclose(full_fail_file);
  } else {
      for (int ii=0;ii<opt.threads;++ii) {
	std::stringstream current_mapfail_file;
	current_mapfail_file << opt.out << "_thread" << ii << "_mapfails.gz";
	if ( remove(current_mapfail_file.str().c_str()) != 0) {
	  perror( "Error deleting mapfail file");
	}
	std::stringstream current_likefail_file;
	current_likefail_file << opt.out << "_thread" << ii << "_likefails.gz";
	if ( remove(current_likefail_file.str().c_str()) != 0) {
	  std::cerr << current_likefail_file.str().c_str() << " "; 
	  perror( "Error deleting likefail file");
	}
      }
  }

  if (opt.likeplot) {

     std::stringstream pre_full_maxlike_file;
     pre_full_maxlike_file << opt.out << ".maxlikes.gz";
     gzFile full_maxlike_file = gzopen(pre_full_maxlike_file.str().c_str(),"wb");
     
     uint64_t maxlikebuffsize = 1ULL<<20;
     char *maxlikebuff = new char[maxlikebuffsize];

    for (int ii=0;ii<opt.threads;++ii) {
      std::stringstream max_likes_file;
      max_likes_file << opt.out << "_thread" << ii << "_maxlikes.gz";
      gzFile infile = gzopen(max_likes_file.str().c_str(),"rb");
      if (!infile) {
	std::cerr << "Unable to open maxlike file for thread " << ii << std::endl;
	exit(1);
      }
      int still_reading = 0;
      while ((still_reading = gzread(infile,maxlikebuff,maxlikebuffsize)) > 0 ) {
	std::string active = std::string(maxlikebuff);
	std::string new_active = active.substr(0,still_reading);
	gzwrite(full_maxlike_file,new_active.c_str(),new_active.size());
      }
      gzclose(infile);	
      if ( remove(max_likes_file.str().c_str()) !=0) {
	perror("Error deleting empty maximum likelihoods file");
      }
    }

    delete[] maxlikebuff;
    gzclose(full_maxlike_file);
  } else {
    for (int ii=0;ii<opt.threads;++ii) {
      std::stringstream max_likes_file;
      max_likes_file << opt.out << "_thread" << ii << "_maxlikes.gz";
      if ( remove(max_likes_file.str().c_str()) !=0) {
	perror("Error deleting empty maximum likelihoods file");
      }
    }
  }
 
	//  like_fails << out_base << "_thread" << thread_num << "_likefails.gz";

  std::cerr << " done" << std::endl;
  logfile_track << " done" << std::endl;

  if (!tc.lenient) {
    std::cerr << "Processed " << pretty_num(numreads) << " reads, "
	      << pretty_num(nummapped) << " reads pseudoaligned" << std::endl;
    logfile_track << "Processed " << pretty_num(numreads) << " reads, "
	    << pretty_num(nummapped) << " reads pseudoaligned" << std::endl;
    nleft = nummapped - like_filter.number_removed;

    if (like_filter.use==true) {
      std::cerr << like_filter.number_removed << " reads removed for failing likelihood filter: z-score below " << like_filter.threshold << ", " << nleft << " reads remain for analysis" << std::endl;
      logfile_track << like_filter.number_removed << " reads removed for failing likelihood filter: z-score below " << like_filter.threshold << ", " << nleft << " reads remain for analysis" << std::endl;
    } else {
      std::cerr << "No likelihood filter, " << nleft << " reads used for analysis" << std::endl;
      logfile_track << "No likelihood filter, " << nleft << " reads used for analysis" << std::endl;
    }
  


} else {
    std::cerr << "Processed " << pretty_num(numreads) << " reads, " <<
      pretty_num(numstrict) << " strictly pseudomapped and " << pretty_num(numresc) << " matched on lenient setting, for a total of " << pretty_num(nummapped) << " reads pseudoaligning" << std::endl;
    logfile_track << "Processed " << pretty_num(numreads) << " reads, " <<
      pretty_num(numstrict) << " strictly pseudomapped and " << pretty_num(numresc) << " matched on lenient setting, for a total of " << pretty_num(nummapped) << " reads pseudoaligning" << std::endl;
    nleft = nummapped - like_filter.number_removed;
    if (like_filter.use==true) {
      std::cerr << like_filter.number_removed << " reads removed for failing likelihood filter: z-score below " << like_filter.threshold << ", " << nleft << " reads remain for analysis" << std::endl;
      logfile_track << like_filter.number_removed << " reads removed for failing likelihood filter: z-score below " << like_filter.threshold << ", " << nleft << " reads remain for analysis" << std::endl;
    } else {
      std::cerr << "No likelihood filter, " << nleft << " reads used for analysis" << std::endl;
      logfile_track << "No likelihood filter, " << nleft << " reads used for analysis" << std::endl;
    }
  }

  // auto diff1 = std::chrono::system_clock::now() - tp1;
  // std::cout << "In Master  " << std::chrono::duration_cast<std::chrono::seconds>(diff1).count() << std::endl;

  
  return nleft;
}

/** -- read processors -- **/

void MasterProcessor::processReads() {
  
  // start worker threads
  std::vector<std::thread> workers;
  for (int i = 0; i < opt.threads; i++) {  
    workers.emplace_back(std::thread(ReadProcessor(index,opt,tc,*this,findex,like_filter,logfile_track,i,earlytax,collapse)));
  }

  // let the workers do their thing
  for (int i = 0; i < opt.threads; i++) {
    workers[i].join(); //wait for them to finish
  }

  //diff = std::chrono::system_clock::now() - tp1;
  //std::cout << "Time to processReads:  " << std::chrono::duration_cast<chrono::seconds>(diff).count() << " sec" << std::endl;

  
}

ReadProcessor::ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, fastaIndex& findex,LKFilter& like_filter,std::stringstream& logfile_track, int thread_num, EarlyTaxonomy& earlytax,bool collapse) :
  paired(opt.single_end), tc(tc), index(index), mp(mp), earlytax(earlytax), findex(findex), like_filter(like_filter), illumina_version(opt.illumina_version), fail(opt.fail), logfile_track(logfile_track), thread_num(thread_num),min_logl(opt.min_logl), out_base(opt.out), collapse(collapse), likeplot(opt.likeplot) {
  //initialize buffer

  bufsize = 1ULL<<22;
  int limit = 1ULL<<22;

  //bufsize = 1ULL<<10;
  //int limit = 1ULL<<10;

  //std::cerr << "Buffer size is " << limit << std::endl;
  
  buffer = new char[bufsize];

  //newEcs.reserve(1000);
  //counts.reserve((int) (tc.counts.size() * 1.25));
  clear();
}

ReadProcessor::~ReadProcessor() {
  //if (buffer) {
    //   delete[] buffer;
    //}
}

void ReadProcessor::operator()() {

  std::stringstream full_outfilename;
  full_outfilename << out_base << "_thread" << thread_num << ".lks";
  gzFile full_outfile = gzopen(full_outfilename.str().c_str(),"wb");

  std::stringstream singles_outfilename;
  singles_outfilename << out_base << "_thread" << thread_num << "_singles.lks";
  gzFile single_outfile = gzopen(singles_outfilename.str().c_str(),"wb");

  std::stringstream map_fails; std::stringstream like_fails;
  map_fails << out_base << "_thread" << thread_num << "_mapfails.gz";
  like_fails << out_base << "_thread" << thread_num << "_likefails.gz";
  gzFile mapfail = gzopen(map_fails.str().c_str(),"wb");
  gzFile likefail = gzopen(like_fails.str().c_str(),"wb");

  std::stringstream max_likes;
  max_likes << out_base << "_thread" << thread_num << "_maxlikes.gz";
  gzFile maxlikes = gzopen(max_likes.str().c_str(),"wb");

  if (!full_outfile || !single_outfile) {
    std::cerr << "Error, unable to open intermediate output files for thread " << thread_num << std::endl;
    exit(1);
  }
  bool keep_going = true;
      
  while (keep_going) {
    // grab the reader lock
    {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR.empty()) {
        // nothing to do
        keep_going = false;
      } else {
        // get new sequences
        mp.SR.fetchSequences(buffer, bufsize, seqs, names, quals);

      }
    } // release the reader lock

    if (keep_going) {
      // process our sequences
 
      thread_output thread_out = processBuffer();

      WriteInt(thread_out.clikes,full_outfile,single_outfile,logfile_track,like_filter,thread_out.c_z_support,likefail,thread_out.read_order,maxlikes);
      for (auto it=thread_out.c_z_support.begin();it!=thread_out.c_z_support.end();++it) {
	auto get = mp.zsupport[thread_num].find(it->first);
	if (get==mp.zsupport[thread_num].end()) {
	  std::pair<unsigned int, bool> nent(it->first,it->second);
	  mp.zsupport[thread_num].insert(nent);
	} else {
	  if (get->second==false && it->second==true) {
	    get->second=true;
	  }
	}

      }
      mp.numreads[thread_num]+=thread_out.numreads;
      mp.numstrict[thread_num]+=thread_out.numstrict;
      mp.numresc[thread_num]+=thread_out.numresc;
      if (fail) {
	gzwrite(mapfail,thread_out.mapfail.c_str(),thread_out.mapfail.size());
      }
    }    
    clear();    
  }
  gzclose(full_outfile);
  gzclose(single_outfile);
  gzclose(mapfail);
  gzclose(likefail);
  gzclose(maxlikes);

  return;
}


std::string ReadProcessor::getSeqString(int& ient) {

  int filenum = findex.values[ient].file;

  std::ifstream fastafile(findex.transfasta[filenum].c_str(), std::ios::in);
  if (!fastafile) {
    std::cerr << "Error, couldn't open fasta file " <<  findex.values[ient].file << std::endl;
    logfile_track << "Error, couldn't open fasta file " <<  findex.values[ient].file << std::endl;
  }

  fastafile.seekg(findex.values[ient].offset, std::ios::beg);
  std::string seqstring;
  getline(fastafile,seqstring);  

  return seqstring;
}

double ReadProcessor::Alignment(const char* &s1,const char* &s2,const char* &q1, const char* &q2,int& ient,const uint64_t& linelength, bool& paired,rt_entry& map_entry,std::string& seqstring) {

  //align the forward read
  StripedSmithWaterman::Aligner aligner;
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment forward_alignment;
  
  aligner.Align(s1, seqstring.c_str(), linelength, filter, &forward_alignment);
  //deal with reverse read
  
  //calculate the likelihoods
  double forward_likelihood = calcLikelihood(forward_alignment.cigar_string,q1,illumina_version);
  
  //deal with reverse read  
  StripedSmithWaterman::Alignment reverse_alignment;

  int s2_size = std::strlen(s2);
  std::string rs2;	  
  for (int jj=0;jj<s2_size;++jj) {
    rs2 = rs2 + twinletter(s2[s2_size-(jj+1)]);	     
  }
    
  aligner.Align(rs2, seqstring.c_str(), linelength, filter, &reverse_alignment);

  double reverse_likelihood = calcLikelihood(reverse_alignment.cigar_string,q2,illumina_version);
  
  double read_likelihood = forward_likelihood + reverse_likelihood;
  auto f1 = map_entry.counts.find(read_likelihood);
  if (f1==map_entry.counts.end()) {
    std::pair<double, int> centry(read_likelihood,1);
    map_entry.counts.insert(centry);
  } else {
    f1->second++;
  }

  std::string seen_string = seqstring.substr(forward_alignment.ref_begin,std::strlen(s1)) + seqstring.substr(reverse_alignment.ref_begin,std::strlen(s2));
  std::pair<std::string, double> nentry(seen_string, read_likelihood);
  map_entry.seen_likes.insert(nentry);
  map_entry.single_starts.push_back(forward_alignment.ref_begin);
  map_entry.paired_starts.push_back(reverse_alignment.ref_begin);
  return read_likelihood;

}

double ReadProcessor::Alignment(const char* &s1,const char* &q1,int& ient,const uint64_t& linelength, bool& paired,rt_entry& map_entry,std::string& seqstring) {

  //align the forward read
  StripedSmithWaterman::Aligner aligner;
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment forward_alignment;
  
  aligner.Align(s1, seqstring.c_str(), linelength, filter, &forward_alignment);
  
  //calculate the likelihoods
  double read_likelihood = calcLikelihood(forward_alignment.cigar_string,q1,illumina_version);
  auto f1 = map_entry.counts.find(read_likelihood);
  if (f1==map_entry.counts.end()) {
    std::pair<double, int> centry(read_likelihood,1);
    map_entry.counts.insert(centry);
  } else {
    f1->second++;
  }
  //std::cout << " " << read_likelihood << std::endl;
  std::string seen_string = seqstring.substr(forward_alignment.ref_begin,std::strlen(s1));
  std::pair<std::string, double> nentry(seen_string, read_likelihood);
  map_entry.seen_likes.insert(nentry);
  map_entry.single_starts.push_back(forward_alignment.ref_begin);
  return read_likelihood;

}

void ReadProcessor::WriteInt(std::vector<std::unordered_map<double, std::vector<unsigned int> > >& clikes,gzFile& full_outfile,gzFile& single_outfile,std::stringstream& logfile_track,LKFilter& like_filter,std::unordered_map<unsigned int, bool>& c_z_support,gzFile& likefail,std::unordered_map<int, std::string>& read_order,gzFile& maxlikes) {

  double thresh_val = like_filter.mean + like_filter.sd*like_filter.threshold;

  int removed = 0;
  std::stringstream full_ss;
  std::stringstream single_ss;
  std::stringstream likefail_ss;
  std::stringstream maxlikes_ss;

  if (like_filter.use==true) {

    for (int i=0;i<clikes.size();++i) {
    
      double max_val = 9;
      for (auto get=clikes[i].begin();get!=clikes[i].end();++get) {
	if (max_val==9) {
	  max_val = get->first;
	} else {
	  if (get->first > max_val) {
	    max_val = get->first;
	  }
	}
      }
      maxlikes_ss << max_val << " ";
      
      if (max_val > thresh_val) {
	if (clikes[i].size() > 1) {
	  for (auto get=clikes[i].begin();get!=clikes[i].end();++get) {	
	    for (int j=0;j< get->second.size();++j) {
	      auto set = c_z_support.find(get->second[j]);
	      if (set==c_z_support.end()) {
		std::cerr << "Error finding " << get->second[j] << " in zero support" << std::endl;
	      } else {
		set->second = true;
	      }	  
	      full_ss << get->first << "|" << get->second[j] << "@";
	    }
	  }	
	  full_ss << "&";
	} else {
	  auto get = clikes[i].begin();
	  if (get->second.size() > 1) {
	    for (int j=0;j< get->second.size();++j) {
	      full_ss << get->first << "|" << get->second[j] << "@";
	      auto set = c_z_support.find(get->second[j]);
	      if (set==c_z_support.end()) {
		std::cerr << "Error finding " << get->second[j] << " in zero support" << std::endl;
	      } else {
		set->second = true;
	  
	      }
	    }
	    full_ss << "&";
	  } else {
	    single_ss << get->second[0] << "&";
	    auto set = c_z_support.find(get->second[0]);
	    if (set==c_z_support.end()) {
	      std::cerr << "Error finding " << get->second[0] << " in zero support" << std::endl;
	    } else {
	      set->second = true;
	    }
	  }
	}
      } else {
	if (fail) {
	  auto read_name = read_order.find(i);
	  if (read_name!=read_order.end()) {
	    likefail_ss << read_name->second << " ";
	  }
	}	
	++removed;
      }
    }
 
    like_filter.number_removed+=removed;  
    gzwrite(full_outfile,full_ss.str().c_str(),full_ss.str().size());
    gzwrite(single_outfile,single_ss.str().c_str(),single_ss.str().size());
    if (fail) {
      gzwrite(likefail,likefail_ss.str().c_str(),likefail_ss.str().size());
    }
    if (likeplot) {
      gzwrite(maxlikes,maxlikes_ss.str().c_str(),maxlikes_ss.str().size());
    }

  } else {

    for (int i=0;i<clikes.size();++i) {

      double max_val = 9;

      if (clikes[i].size() > 1) {
	for (auto get=clikes[i].begin();get!=clikes[i].end();++get) {
	  if (max_val==9) {
	    max_val = get->first;
	  } else {
	    if (get->first > max_val) {
	      max_val = get->first;
	    }
	  }
	  for (int j=0;j< get->second.size();++j) {
	    auto set = c_z_support.find(get->second[j]);
	    if (set==c_z_support.end()) {
	      std::cerr << "Error finding " << get->second[j] << " in zero support" << std::endl;
	    } else {
	      set->second = true;
	    }	  
	    full_ss << get->first << "|" << get->second[j] << "@";
	  }
	}	
	full_ss << "&";
      } else {
	auto get = clikes[i].begin();
	max_val = get->first;
	if (get->second.size() > 1) {
	  for (int j=0;j< get->second.size();++j) {
	    full_ss << get->first << "|" << get->second[j] << "@";
	    auto set = c_z_support.find(get->second[j]);
	    if (set==c_z_support.end()) {
	      std::cerr << "Error finding " << get->second[j] << " in zero support" << std::endl;
	    } else {
	      set->second = true;
	  
	    }
	  }
	  full_ss << "&";
	} else {
	  single_ss << get->second[0] << "&";
	  auto set = c_z_support.find(get->second[0]);
	  if (set==c_z_support.end()) {
	    std::cerr << "Error finding " << get->second[0] << " in zero support" << std::endl;
	  } else {
	    set->second = true;
	  }
	}
      }
      maxlikes_ss << max_val << " ";
     
    }
    gzwrite(full_outfile,full_ss.str().c_str(),full_ss.str().size());
    gzwrite(single_outfile,single_ss.str().c_str(),single_ss.str().size());
    if (likeplot) {
      gzwrite(maxlikes,maxlikes_ss.str().c_str(),maxlikes_ss.str().size());
    }
  }
}

thread_output ReadProcessor::processBuffer() {

  // set up thread variables
  std::vector<std::pair<KmerEntry,int>> v1, v2;
  std::vector<int> vtmp;
  std::vector<int> u;
  int inlike_row = 0;
  std::vector<std::unordered_map<double, std::vector<unsigned int> > > current_likelihoods;
  uint32_t numresc = 0;
  uint32_t numstrict = 0;
  std::unordered_map<unsigned int,bool > in_z;
  std::stringstream mapfail_ss;
  std::unordered_map<int, std::string> read_order;
  int current_mapped = 0;

  u.reserve(1000);
  v1.reserve(1000);
  v2.reserve(1000);
  vtmp.reserve(1000);

  const char* s1 = 0;
  const char* s2 = 0;
  int l1,l2;

  // std::chrono::system_clock::time_point tp1;
  // std::chrono::system_clock::time_point tp2;
  // tp1 = std::chrono::system_clock::now();
  // auto diff1 = std::chrono::system_clock::now() - tp1;
  
  // actually process the sequences
  for (int i = 0; i < seqs.size(); i++) {
    
    s1 = seqs[i].first;
    l1 = seqs[i].second;
    if (paired) {
      i++;
      s2 = seqs[i].first;
      l2 = seqs[i].second;
    }

    v1.clear();
    v2.clear();
    u.clear();

  
    // process read
    index.match(s1,l1, v1);
    if (paired) {
      index.match(s2,l2, v2);
    }

    // collect the target information
    
    int ec = -1;
    int r = tc.intersectKmers(v1, v2, !paired,u);
    if (u.empty()) {
      if (tc.lenient) {
	int ret = tc.altIntersectKmers(v1,v2,!paired,u);
	if (u.empty()) {
	  if (fail) {
	    mapfail_ss << std::string(names[i].first,strlen(names[i].first)) << ",";
	  }
	  continue;
	} else {
	  ++numresc;
	}
      } else {
	if (fail) {
	  mapfail_ss << std::string(names[i].first,strlen(names[i].first)) << ",";
	}
    	continue;
      }
    } else {
      ++numstrict;
    }
 

    if (collapse==true) {

      read_tmap tmap;
      
      for (auto tr : u) {

	std::unordered_map<std::string,int>::const_iterator got = findex.ref_names.find(index.target_names_[tr].c_str());
	if ( got == findex.ref_names.end() ) {
	  std::cerr << index.target_names_[tr]  << "  not found in fasta reference, is index up to date?\n";
	  logfile_track << index.target_names_[tr]  << "  not found in fasta reference, is index up to date?\n";
	} else {
	  int ient = got->second;
	  auto tlook = earlytax.label_map.find(ient);
	  if (tlook==earlytax.label_map.end()) {
	    unsigned int ladd = earlytax.label_map.size() + 1;
	    std::pair<int, unsigned int> newlab(ient,ladd);
	    earlytax.label_map.insert(newlab);
	    std::pair<unsigned int, unsigned int> newsize(ladd,1);
	    earlytax.size_map.insert(newsize);
	    std::pair<std::string, unsigned int> newname(got->first,ladd);
	    earlytax.name_map.insert(newname);
	    tlook = earlytax.label_map.find(ient);
	  }
	  //std::cout << got->first << " " << tlook->second << std::endl;
	  auto rget = tmap.find(tlook->second);
	  
	  if (rget==tmap.end()) {
	    
	    rt_entry new_entry;
	    //std::cout << "New entry for " << tlook->second << " values ";

	    //tp1 = std::chrono::system_clock::now();
	    std::string seqstring = getSeqString(ient);
	    //tp2 = std::chrono::system_clock::now();
	    //diff1+=tp2-tp1;
	    //std::cout << "<" << got->first << std::endl << seqstring << std::endl;
	    if (paired) {
	      Alignment(s1,s2,quals[i-1].first,quals[i].first,ient, findex.values[ient].linelength,paired,new_entry,seqstring);
	    } else {
	      Alignment(s1,quals[i].first,ient, findex.values[ient].linelength,paired,new_entry,seqstring);
	    }
	    std::pair<unsigned int, rt_entry> tmap_entry(tlook->second,new_entry);
	    tmap.insert(tmap_entry);
	    
	   	    
	  } else {
    
	    //	    tp1 = std::chrono::system_clock::now();
	    std::string seqstring = getSeqString(ient);
	    //	    tp2 = std::chrono::system_clock::now();
	    //diff1+=tp2-tp1;
	    bool found = false;
	    //std::cout << "For taxa " << tlook->second;
	    for (int ii=0;ii<rget->second.single_starts.size();++ii) {
	      bool check = false;
	      if (paired) {
		if ((rget->second.single_starts[ii] + std::strlen(s1)) < seqstring.length() && (rget->second.paired_starts[ii] + std::strlen(s2)) < seqstring.length()) {
		  check = true;
		} 
	      } else {
		if ((rget->second.single_starts[ii] + std::strlen(s1)) < seqstring.length()) {
		  check = true;
		}
	      }         
	      if (check==true) {
		std::string check_string;
		if (paired) {
		  check_string = seqstring.substr(rget->second.single_starts[ii],std::strlen(s1)) + seqstring.substr(rget->second.paired_starts[ii],std::strlen(s2));
		} else {
		  check_string = seqstring.substr(rget->second.single_starts[ii],std::strlen(s1));
		}
		auto sget = rget->second.seen_likes.find(check_string);
		if (sget!=rget->second.seen_likes.end()) {
		  auto cget = rget->second.counts.find(sget->second);
		  //std::cout << " " << sget->second << std::endl;
		  if (cget!=rget->second.counts.end()) {
		    cget->second++;
		    found = true;
		    break;
		  } else {
		    std::cerr << "Error with read map for read " << s1 << " and taxa " << tlook->second << std::endl;
		    exit(1);
		  }
		}
	      }
	    }
	    
	    if (found==false) {
	      //std::cout << " need to align ";
	      if (paired) {
		Alignment(s1,s2,quals[i-1].first,quals[i].first,ient, findex.values[ient].linelength,paired,rget->second,seqstring);
	      } else {
		Alignment(s1,quals[i].first,ient, findex.values[ient].linelength,paired,rget->second,seqstring);
	      }
	      
	    }
	  }
	}
      }

      std::unordered_map<double, std::vector<unsigned int>> read_likes = interpretTmap(tmap);
      for (auto it = read_likes.begin();it!=read_likes.end();++it) {
      	std::vector<unsigned int> seenvec = it->second;
      	for (int ii=0;ii<seenvec.size();++ii) {
      	  auto finder = in_z.find(seenvec[ii]);
      	  if (finder==in_z.end()) {
      	    std::pair<unsigned int, bool> news(seenvec[ii],false);
      	    in_z.insert(news);
      	  }
      	}
      }
      current_likelihoods.push_back(read_likes);


    } else {

      /////////////////////////////

      //std::cerr << "READ " << i << " size " << u.size() << std::endl;
      
      read_tmap tmap;
      std::unordered_map<double, std::vector<unsigned int>> read_likes;

	for (auto tr : u) {

	  std::unordered_map<std::string,int>::const_iterator got = findex.ref_names.find(index.target_names_[tr].c_str());
	  if ( got == findex.ref_names.end() ) {
	    std::cerr << index.target_names_[tr]  << "  not found in fasta reference, is index up to date?\n";
	    logfile_track << index.target_names_[tr]  << "  not found in fasta reference, is index up to date?\n";
	  } else {
	    int ient = got->second;
	    unsigned int temp_ient = (unsigned int) ient;
	    auto tlook = earlytax.label_map.find(ient);
	    if (tlook==earlytax.label_map.end()) {
	      unsigned int ladd = earlytax.label_map.size() + 1;
	      std::pair<int, unsigned int> newlab(ient,ladd);
	      earlytax.label_map.insert(newlab);
	      std::pair<unsigned int, unsigned int> newsize(ladd,1);
	      earlytax.size_map.insert(newsize);
	      std::pair<std::string, unsigned int> newname(got->first,ladd);
	      earlytax.name_map.insert(newname);
	      tlook = earlytax.label_map.find(ient);
	    }
	    auto rget = tmap.find(tlook->second);
	  
	    if (rget==tmap.end()) {
	    
	      rt_entry new_entry;
	      std::string seqstring = getSeqString(ient);

	      if (paired) {
		double likelihood = Alignment(s1,s2,quals[i-1].first,quals[i].first,ient, findex.values[ient].linelength,paired,new_entry,seqstring);
		auto growl = read_likes.find(likelihood);
		if (growl==read_likes.end()) {
		  std::vector<unsigned int> temp_vec;
		  temp_vec.push_back(temp_ient);
		  std::pair<double, std::vector<unsigned int > > clike_entry(likelihood,temp_vec);
		  read_likes.insert(clike_entry);
		} else {
		  growl->second.push_back(temp_ient);
		}
	      } else {
		double likelihood = Alignment(s1,quals[i].first,ient, findex.values[ient].linelength,paired,new_entry,seqstring);
		auto growl = read_likes.find(likelihood);
		if (growl==read_likes.end()) {
		  std::vector<unsigned int> temp_vec;
		  temp_vec.push_back(temp_ient);
		  std::pair<double, std::vector<unsigned int > > clike_entry(likelihood,temp_vec);
		  read_likes.insert(clike_entry);
		} else {
		  growl->second.push_back(temp_ient);
		}
	      }
	      std::pair<unsigned int, rt_entry> tmap_entry(tlook->second,new_entry);
	      tmap.insert(tmap_entry);
	    	   	    
	    } else {
    
	      std::string seqstring = getSeqString(ient);
	      bool found = false;
	      for (int ii=0;ii<rget->second.single_starts.size();++ii) {
		bool check = false;
		if (paired) {
		  if ((rget->second.single_starts[ii] + std::strlen(s1)) < seqstring.length() && (rget->second.paired_starts[ii] + std::strlen(s2)) < seqstring.length()) {
		    check = true;
		  } 
		} else {
		  if ((rget->second.single_starts[ii] + std::strlen(s1)) < seqstring.length()) {
		    check = true;
		  }
		}
		if (check==true) {
		  std::string check_string;
		  if (paired) {
		    		    		    
		    check_string = seqstring.substr(rget->second.single_starts[ii],std::strlen(s1)) + seqstring.substr(rget->second.paired_starts[ii],std::strlen(s2));
		    
		  } else {
		    check_string = seqstring.substr(rget->second.single_starts[ii],std::strlen(s1));
		  }
		  auto sget = rget->second.seen_likes.find(check_string);
		  if (sget!=rget->second.seen_likes.end()) {


		    auto growl = read_likes.find(sget->second);
		    if (growl==read_likes.end()) {
		      std::cerr << "Error with seen likelihood not appearing in hash: " << rget->first;
		      exit(1);
		    }
		    growl->second.push_back(temp_ient);

		    auto cget = rget->second.counts.find(sget->second);
		    if (cget!=rget->second.counts.end()) {
		      cget->second++;
		      found = true;
		      break;
		    } else {
		      std::cerr << "Error with read map for read " << s1 << " and taxa " << tlook->second << std::endl;
		      exit(1);
		    }
		  }
		}
	      }
	    
	      if (found==false) {
		if (paired) {
		  double likelihood = Alignment(s1,s2,quals[i-1].first,quals[i].first,ient, findex.values[ient].linelength,paired,rget->second,seqstring);
		  auto growl = read_likes.find(likelihood);
		  if (growl==read_likes.end()) {
		    std::vector<unsigned int> temp_vec;
		    temp_vec.push_back(temp_ient);
		    std::pair<double, std::vector<unsigned int > > clike_entry(likelihood,temp_vec);
		    read_likes.insert(clike_entry);
		  } else {
		    growl->second.push_back(temp_ient);
		  }
		} else {
		  double likelihood = Alignment(s1,quals[i].first,ient, findex.values[ient].linelength,paired,rget->second,seqstring);
		  auto growl = read_likes.find(likelihood);
		  if (growl==read_likes.end()) {
		    std::vector<unsigned int> temp_vec;
		    temp_vec.push_back(temp_ient);
		    std::pair<double, std::vector<unsigned int > > clike_entry(likelihood,temp_vec);
		    read_likes.insert(clike_entry);
		  } else {
		    growl->second.push_back(temp_ient);
		  }
		}
	      }
	    }
	  }
	}
    
      //std::unordered_map<double, std::vector<unsigned int>> read_likes = interpretTmap(tmap);
      for (auto it = read_likes.begin();it!=read_likes.end();++it) {
      	std::vector<unsigned int> seenvec = it->second;
      	for (int ii=0;ii<seenvec.size();++ii) {
      	  auto finder = in_z.find(seenvec[ii]);
      	  if (finder==in_z.end()) {
      	    std::pair<unsigned int, bool> news(seenvec[ii],false);
      	    in_z.insert(news);
      	  }
      	}
      }
      current_likelihoods.push_back(read_likes);

      //////////////////////////////
    }
    
    if (fail) {
      std::pair<int, std::string> trackorder(current_mapped,std::string(names[i].first,strlen(names[i].first)));
      read_order.insert(trackorder);
      current_mapped++;
    }

  }
    
  thread_output current_output;
  current_output.c_z_support = in_z;
  current_output.numresc = numresc;
  current_output.numstrict = numstrict;
  if (!paired) {
    current_output.numreads = seqs.size();
  } else {
    current_output.numreads = seqs.size()/2;
  }
  current_output.clikes = current_likelihoods;
  current_output.mapfail = mapfail_ss.str();
  current_output.read_order = read_order;
  //current_output.c_hf_qualityscores = current_qcs;
  //current_output.c_z_support = current_z_support;
  return current_output;  
  
}

std::unordered_map<double, std::vector<unsigned int>> ReadProcessor::interpretTmap(read_tmap& tmap) {

  std::unordered_map<double, std::vector<unsigned int>> in_read_like;

    for (auto it=tmap.begin();it!=tmap.end();++it) {
      double like = 0;
      auto ts = earlytax.size_map.find(it->first);
      if (ts==earlytax.size_map.end()) {
	std::cerr << "Error with earlytaxa map\n";
	exit(1);
      }
      double tot_num = (double) ts->second;
      double remaining = tot_num;
      //std::cout << "Taxa: " << it->first << " ";
      for (auto ic=it->second.counts.begin();ic!=it->second.counts.end();++ic) {
	double num_seen = (double) ic->second;
	like+=num_seen * exp(ic->first);
	remaining-=num_seen;
	//std::cout << " " << like << " " << ic->first << "(" << num_seen << ") ";
      }
      like+=remaining*exp(min_logl);
      like = log(like/tot_num);
      //std::cout << "F: " << like << std::endl;
      auto sset = in_read_like.find(like);
      if (sset==in_read_like.end()) {
	std::vector<unsigned int> intype;
	intype.push_back(it->first);
	std::pair<double, std::vector<unsigned int> > extentry(like,intype);
	in_read_like.insert(extentry);
      } else {
	sset->second.push_back(it->first);
      }
    }

  return in_read_like;
}

void ReadProcessor::clear() {
  memset(buffer,0,bufsize);
  //newEcs.clear();
  //counts.clear();
  //counts.resize(tc.counts.size(),0);
}

/** -- sequence reader -- **/
SequenceReader::~SequenceReader() {
  if (fp1) {
    gzclose(fp1);
  }
  if (paired && fp2) {
    gzclose(fp2);
  }

  kseq_destroy(seq1);
  if (paired) {
    kseq_destroy(seq2);
  }
}

// returns true if there is more left to read from the files
bool SequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
  std::vector<std::pair<const char *, int> > &names,
  std::vector<std::pair<const char *, int> > &quals) {
  
  seqs.clear();
  names.clear();
  quals.clear();

  int bufpos = 0;
  int pad = (paired) ? 2 : 1;
  while (true) {
    if (!state) { // should we open a file
      if (current_file >= files.size()) {
        // nothing left
        return false;
      } else {
        // close the current file
        if(fp1) {
          gzclose(fp1);
        }
        if (paired && fp2) {
          gzclose(fp2);
        }
        // open the next one
        fp1 = gzopen(files[current_file].c_str(),"r");
        seq1 = kseq_init(fp1);
        l1 = kseq_read(seq1);
        state = true;
        if (paired) {
          current_file++;
          fp2 = gzopen(files[current_file].c_str(),"r");
          seq2 = kseq_init(fp2);
          l2 = kseq_read(seq2);
        }
      }
    }
    // the file is open and we have read into seq1 and seq2

    if (l1 > 0 && (!paired || l2 > 0)) {
      int bufadd = l1 + l2 + pad;
      // fits into the buffer
      //if (full) {
        nl1 = seq1->name.l;
        if (paired) {
          nl2 = seq2->name.l;
        }
        bufadd += (l1+l2) + pad + (nl1+nl2)+ pad;
	//}
      if (bufpos+bufadd< limit) {
        char *p1 = buf+bufpos;
        memcpy(p1, seq1->seq.s, l1+1);
        bufpos += l1+1;
        seqs.emplace_back(p1,l1);
        //if (full) {
          p1 = buf+bufpos;
          memcpy(p1, seq1->qual.s,l1+1);
          bufpos += l1+1;
          quals.emplace_back(p1,l1);
          p1 = buf+bufpos;
          memcpy(p1, seq1->name.s,nl1+1);
          bufpos += nl1+1;
          names.emplace_back(p1,nl1);
	  //}

        if (paired) {
          char *p2 = buf+bufpos;
          memcpy(p2, seq2->seq.s,l2+1);
          bufpos += l2+1;
          seqs.emplace_back(p2,l2);
	  //  if (full) {
            p2 = buf+bufpos;
            memcpy(p2,seq2->qual.s,l2+1);
            bufpos += l2+1;
            quals.emplace_back(p2,l2);
            p2 = buf + bufpos;
            memcpy(p2,seq2->name.s,nl2+1);
            bufpos += nl2+1;
            names.emplace_back(p2,nl2);
	    //   }
        }
      } else {
        return true; // readit next time
      }

      // read for the next one
      l1 = kseq_read(seq1);
      if (paired) {
        l2 = kseq_read(seq2);
      }
    } else {
      current_file++; // move to next file
      state = false; // haven't opened file yet
    }
  }
}

bool SequenceReader::empty() {
  return (!state && current_file >= files.size());
}

