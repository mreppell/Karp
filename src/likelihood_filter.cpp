
#include "likelihood_filter.hpp"
#include "ProcessReads.h"


void LKFilter::justGetReadNum(const ProgramOptions& opt) {

  SequenceReader SR(opt);
  size_t bufsize = 1ULL<<20;
  char* buffer = new char[bufsize];
  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  bool keep_going = true;
  while (keep_going) {
    if (SR.empty()) {
      // nothing to do
      keep_going = false;
    } else {
      // get new sequences
      SR.fetchSequences(buffer, bufsize, seqs, names, quals);
      numreads+=seqs.size();
    }

    memset(buffer,0,bufsize);
  } 

  seqs.clear();
  quals.clear();
  names.clear();
  delete[] buffer;

}

void LKFilter::GetParameterValues(const ProgramOptions& opt) {


  SequenceReader SR(opt);
  ArrayXXi hf_qualityscores(42,1);

  size_t bufsize = 1ULL<<20;
  char* buffer = new char[bufsize];
  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  bool keep_going = true;

  while (keep_going) {
    if (SR.empty()) {
      // nothing to do
      keep_going = false;
    } else {
      // get new sequences
      SR.fetchSequences(buffer, bufsize, seqs, names, quals);
      numreads+=seqs.size();
      lk_processBuffer(hf_qualityscores,quals);
    }

   memset(buffer,0,bufsize);
  } 

  seqs.clear();
  quals.clear();
  names.clear();
  delete[] buffer;

  std::vector<double> results = E_logl(hf_qualityscores);
  mean = results[0];
  sd = sqrt(results[1]);

}

void LKFilter::lk_processBuffer(ArrayXXi& hf_qualityscores,std::vector<std::pair<const char*, int>>& quals) {


  if (!paired) {

    for (int read = 0;read < quals.size(); ++read) {
    
      if (hf_qualityscores.cols()==1) {
	hf_qualityscores.resize(42,strlen(quals[read].first));
	hf_qualityscores.setZero();
      }

      if (strlen(quals[read].first) > hf_qualityscores.cols()) {
	
	ArrayXXi new_hf_qualityscores(42,strlen(quals[read].first));
	new_hf_qualityscores.setZero();
	for (int jj=0;jj<hf_qualityscores.cols();++jj) {
	  for (int ii=0;ii<hf_qualityscores.rows();++ii) {
	    new_hf_qualityscores(ii,jj) = hf_qualityscores(ii,jj);
	  }
	}

	hf_qualityscores = new_hf_qualityscores;
	
      }
      
      updateFilter(hf_qualityscores,quals[read].first,NULL,paired,illumina_version);      
    }
  } else {

    for (int read = 1;read < quals.size(); read+=2) {

      if (hf_qualityscores.cols()==1) {
	hf_qualityscores.resize(42,strlen(quals[read].first)+strlen(quals[read-1].first));
	hf_qualityscores.setZero();
      }

      if (strlen(quals[read].first) + strlen(quals[read-1].first) > hf_qualityscores.cols()) {

	ArrayXXi new_hf_qualityscores(42,strlen(quals[read].first) + strlen(quals[read-1].first));
	new_hf_qualityscores.setZero();
	for (int jj=0;jj<hf_qualityscores.cols();++jj) {
	  for (int ii=0;ii<hf_qualityscores.rows();++ii) {
	    new_hf_qualityscores(ii,jj) = hf_qualityscores(ii,jj);
	  }
	}

	hf_qualityscores = new_hf_qualityscores;
      }

      updateFilter(hf_qualityscores,quals[read-1].first,quals[read].first,paired,illumina_version);
    }
  }
}

//Creates matrix of possible base quality scores by read length, used by harp likelihood filter
void LKFilter::updateFilter(ArrayXXi& hf_qualityscores,const char* quals1,const char* quals2,bool paired,const double& illumina_version) {

  int mod = (int) illumina_version;

  if (!paired) {
    int pos = 0;

    for (auto it = 0;it<strlen(quals1);++it) {
      int row = quals1[it] - mod;
      
      if (row < 0) { std::cerr << "\nError, incorrect illumina quality scores " << quals1[it] << " not consistent with ASCII offset " << mod << std::endl;
	exit(1);
      }      
      if (row > hf_qualityscores.rows() -1) {
	hf_qualityscores.conservativeResize(row + 1,Eigen::NoChange);
      }
      hf_qualityscores(row,pos)++;
      ++pos;
    }
  } else {
    int pos = 0;
    for (auto it = 0;it<strlen(quals1);++it) {
      int row = quals1[it] - mod;
      if (row < 0) { std::cerr << "\nError, incorrect illumina quality scores " << quals1[it] << " not consistent with ASCII offset " << mod << std::endl;
	exit(1);
      }
      if (row > hf_qualityscores.rows() -1) {
	hf_qualityscores.conservativeResize(row + 1,Eigen::NoChange);
      }
      hf_qualityscores(row,pos)++;
      ++pos;
    }
    for (auto it = 0;it<strlen(quals2);++it) {
      int row = quals2[it] - mod;
      if (row < 0) { std::cerr << "\nError, incorrect illumina quality scores " << quals1[it] << " not consistent with ASCII offset " << mod << std::endl;
	exit(1);
      }
      if (row > hf_qualityscores.rows() -1) {
	hf_qualityscores.conservativeResize(row + 1,Eigen::NoChange);
      }
      hf_qualityscores(row,pos)++;
      ++pos;
    }
  }

}

//Calculate expected distribution of read likelihoods for filter

std::vector<double> LKFilter::E_logl_per_site_conditional(size_t q) {

  std::vector<double> result;result.push_back(0);result.push_back(0);

  if (q>0) {

    double qq = (double) q;
    double e = 0.75;
    if (q!=2) {
      e = pow(10, -qq/10);
    }

    result[0] = (1-e)*log(1-e) + e*log(e/3);
    result[1] = (1-e)*log(1-e)*log(1-e) + e*log(e/3)*log(e/3);
  }
  
  return result;
}

std::vector<double> LKFilter::E_logl_per_site(size_t position,ArrayXXi& hf_qualityscores) {

  std::vector<double> result;result.push_back(0);result.push_back(0);

  ArrayXd pos_probs(hf_qualityscores.rows());
  pos_probs.setZero();
  double denom = (double) hf_qualityscores.col(position).sum();
  for (size_t i=0;i<hf_qualityscores.rows();++i) {
    pos_probs(i) = ((double) hf_qualityscores(i,position))/denom;
  }
  for (size_t q=0; q<hf_qualityscores.rows();++q) {
    std::vector<double> in_res = E_logl_per_site_conditional(q);
    result[0]+= in_res[0]*pos_probs(q);
    result[1]+= in_res[1]*pos_probs(q);
     
  }
  return result;

}

std::vector<double> LKFilter::E_logl(ArrayXXi& hf_qualityscores) {
  std::vector<double> result;
  result.push_back(0);result.push_back(0);

  for (size_t position=0; position < hf_qualityscores.cols();++position) {
    std::vector<double> in_res = E_logl_per_site(position,hf_qualityscores);
    result[0]+=in_res[0];
    result[1]+=in_res[1] - (in_res[0]*in_res[0]);
  }
  return result;
}
