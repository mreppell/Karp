#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include<tclap/CmdLine.h>
#include <thread>
#include <time.h>
#include <ctime>
#include <cstdio>
#include <zlib.h>
#include<karpeigen/Eigen/Core>
#include<karpeigen/Eigen/Dense>

#include<chrono>

#include "ssw_cpp.h"
#include "common_k.h"
#include "fastaIndex.h"
#include "ProcessReads.h"
#include "KmerIndex.h"
#include "Kmer.hpp"
#include "MinCollector.h"
#include "karp_like_multi.hpp"
#include "taxonomy.hpp"
#include "likelihood_filter.hpp"
#include "tabulate.hpp"

#define ERROR_STR "Error:"
using namespace std;

std::string asString (const std::chrono::system_clock::time_point& tp)
{
     // convert to system time:
     std::time_t t = std::chrono::system_clock::to_time_t(tp);
     std::string ts = std::ctime(&t);// convert to calendar time
     ts.resize(ts.size()-1);         // skip trailing newline
     return ts;
}

void printConfig(std::ostream& logfile, ProgramOptions& opt) {

  logfile << "Karp " << opt.version << std::endl;
  logfile << "###############Configuration###############\n";
  logfile << "Fastq files\n";
  for (int ii = 0;ii < opt.files.size();++ii) {
    logfile << "\t" << opt.files[ii];
    if (opt.single_end) {
      if (opt.file_direction[ii]) { 
	logfile << " F";
      } else {
	logfile << " R";
      }
    }
    logfile << std::endl;
  }
  logfile << "Reference fasta files\n";
  for (int ii = 0;ii < opt.transfasta.size();++ii) {
    logfile << "\t" << opt.transfasta[ii] << std::endl;
  }
  logfile << "Reference taxonomy files\n";
  for (int ii =0;ii < opt.taxonomy_files.size();++ii) {
    logfile << "\t" << opt.taxonomy_files[ii] << std::endl;
  }
  logfile << "Output files\n\t" << opt.out << ".freqs\t" << opt.out << ".log\n";
  logfile << "Pseudomapping\n";
  logfile << "\tKmer length: " << opt.k << std::endl;
  logfile << "\tIndex file: " << opt.index << std::endl;
  logfile << "EM algorithm\n";
  logfile << "\tBase errors coded according to Illumina version " << opt.illumina_version << std::endl;
  logfile << "\tMin frequency cutoff: " << opt.minimum_frequency_cutoff << std::endl;
  //  logfile << "\tEM restarts: " << opt.em_restarts << std::endl;
  logfile << "\tMax EM iterations: " << opt.max_em_iterations << std::endl;
  logfile << "\tEM Convergence Criteria: " << opt.em_converge << std::endl;
  if (opt.collapse) {
    logfile << "Running in Collapse mode" << std::endl;
  }
  logfile << "Using " << opt.threads << " thread(s)\n\n";
  logfile << "####################Run Log#####################\n";
}

bool checkIndexOptions(ProgramOptions& opt) {

  bool ret = true;

  if (opt.k<=2 || opt.k>= Kmer::MAX_K) {
    cerr << "Error: invalid kmer length " << opt.k << ", min==3 max==" << (Kmer::MAX_K -1) << endl;
    ret=false;
  }

  if (opt.k % 2 == 0) {
    cerr << "Error: k needs to be an odd number" << endl;
    ret = false;
  }

  for (auto& fasta : opt.transfasta) {
    struct stat stFileInfo;
    auto intStat = stat(fasta.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << "Error: Fasta file " << fasta << " not found" << endl;
      ret = false;
    }
  }
  
  if (opt.index=="null") {
    cerr << "Error: no output index file specified" << endl;
    ret = false;
  }
  return ret;
}
    	
int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("For more detailed parameter descriptions please see the users guide", ' ', "salted_caramel_v2.0.1");
    TCLAP::ValueArg<std::string> coMmand("c","command","Command to run, one of \"index,\" \"quantify,\" or \"tabulate\"",true,"null","string");
    cmd.add(coMmand);
    TCLAP::SwitchArg singleStrand("","paired","If classifying separate forward and reverse fastq files from paired-end sequencing, include this flag so Karp knows to expect -q files",false);
    cmd.add(singleStrand);
    TCLAP::ValueArg<std::string> refFasta("r","ref","Comma separated list of reference fasta database files",false,"null","string");
    cmd.add(refFasta);
    TCLAP::ValueArg<std::string> tabSamples("","samples","File with list of sample output files to tabulate with \"tabulate\" command. File should have two columns, with each line a tab or comma separated pair of sample name and corresponding results file",false,"null","string");
    cmd.add(tabSamples);
    TCLAP::ValueArg<std::string> fastQs("f","forward","Comma separated list of fastq files containing reads to quantify, if single end enter fastq filenames with this command\nIf paired end enter only forward direction fastq files with this command and reverse direction files with -q command",false,"null","string");
    cmd.add(fastQs);
    TCLAP::ValueArg<std::string> revQs("q","reverse","If paired-end reads are being classified use this command to enter fastq files that have reverse orientation, should be in same order as partner forward reads entered with -f command. Comma separate multiple files",false,"null","string");
    cmd.add(revQs);
    TCLAP::ValueArg<std::string> kIndex("i","index","De Bruijn index input/output file. If command is \"index\" this is output filename, if \"quantify\" it is an input filename",false,"null","string");
    cmd.add(kIndex);
    TCLAP::ValueArg<int> kmerLength("k","kmer","Length of kmers to use for index",false,31,"int");
    cmd.add(kmerLength);
    TCLAP::ValueArg<int> threadCount("","threads","Number of threads to use",false,1,"int");
    cmd.add(threadCount);
    TCLAP::ValueArg<double> illVersion("","phred","Phred offset for base quality scores, default is Phred+33 corresponding to Illumina 1.8+. Other recognized option is \"64\" for Phred+64",false,33,"double");
    cmd.add(illVersion);
    TCLAP::ValueArg<int> maxIterations("","max_em_it","Maximum number of iterations of EM algorithm before declaring failure to converge [default=250]",false,1000,"int");
    cmd.add(maxIterations);
    TCLAP::ValueArg<double> minRefFreq("","min_freq","Minimum frequency observable for reference haplotypes, values below this are rounded to zero [default=1/Number of Reads]",false,-99,"double");  
    cmd.add(minRefFreq);
    TCLAP::ValueArg<std::string> taxFile("t","tax","Taxonomy file(s) with entries matching reference fasta files. Comma seperated if multiple files. See example folder for format examples.",false,"null","string");
    cmd.add(taxFile);     
    TCLAP::ValueArg<std::string> outFile("o","out","Base name for output files",false,"karp","string");
    cmd.add(outFile);
    TCLAP::SwitchArg collapseTax("","collapse","Collapse taxonomic labels before estimating frequencies [default=FALSE]",false);
    cmd.add(collapseTax);
    TCLAP::SwitchArg readFail("","fail","Output list of read names that fail pseudomapping or likelihood filter [default=FALSE]",false);
    cmd.add(readFail);
    TCLAP::SwitchArg harpFilter("","no_harp_filter","Filter reads using expected distribution of likelihoods? [default=TRUE]",true);
    cmd.add(harpFilter);
    TCLAP::ValueArg<std::string> EMconverg("","em_tolerance","When squared sum of allele frequency changes is lower than this value between algorithm iterations, the EM has converged [default=1e-12]",false,"1e-1","string");
    cmd.add(EMconverg);
    TCLAP::SwitchArg lenientKmers("","strict","Switch arg, use the strict intersection of equivalence classes to declare a pseudomapping match, same condition as Kallisto. If not entered on command line, default behavior is after pseudomapping with strict intersection requirements, reads that fail to map are checked and aligned against references that match multiple but not all kmers. Increases speed but decreases accuracy.",true);
    cmd.add(lenientKmers);
    TCLAP::ValueArg<double> harpFilterThresh("","like_thresh","Likelihood based cutoff. Filter uses z-scores for read likelihoods and filters reads with maximum z-scores below this threshold [default=-2]",false,-2.0,"double");
    cmd.add(harpFilterThresh);
    TCLAP::SwitchArg likePlot("","max_like_out","Output file containing maximum likelihood value for each read, for diagnosing whether to use likelihood filter option",false);
    cmd.add(likePlot);

    cmd.parse(argc,argv);

    ProgramOptions opt;
    std::string command = coMmand.getValue();
    opt.single_end = singleStrand.getValue();
    opt.index = kIndex.getValue();
    opt.k = kmerLength.getValue();
    opt.threads = threadCount.getValue();
    opt.illumina_version = illVersion.getValue();
    opt.minimum_frequency_cutoff = minRefFreq.getValue();
    opt.max_em_iterations = maxIterations.getValue();
    //opt.em_restarts = emRestarts.getValue();
    //opt.format = outFormat.getValue();
    opt.out = outFile.getValue();
    opt.harp_filter = harpFilter.getValue();
    //opt.buffer_size = bufferSize.getValue();
    opt.multi_match = lenientKmers.getValue();
    opt.harp_filter_thresh = harpFilterThresh.getValue();
    opt.collapse = collapseTax.getValue();
    opt.fail = readFail.getValue();
    opt.likeplot= likePlot.getValue();
    //      opt.it_skip = initSkip.getValue();
    std::string all_samples = tabSamples.getValue();
    
    std::string version = "release v1.0";
    opt.version = version;
    
    std::string temp_tol = EMconverg.getValue();
    if (temp_tol.compare("1e-1")==0) {
     	opt.em_converge = 0.000000000001;
    } else {
      opt.em_converge = atof(temp_tol.c_str());
    }


    std::string tax_files = taxFile.getValue();
    std::string rev_fastqs = revQs.getValue();
    std::string ref = refFasta.getValue();
    std::string forward_fastqs = fastQs.getValue();
    
     std::string delim = ",";

     if ((command=="index") || (command=="quantify")) {
       if (ref.compare("null")==0) {
	 std::cerr << "Error, missing reference fasta file, enter with -r\n";
	 exit(1);
       } else {
	 split(ref,delim,opt.transfasta);
       }
     }
     split(tax_files,delim,opt.taxonomy_files);
       

     if (opt.single_end && command=="quantify") {

	std::vector<std::string> forfiles;
	std::vector<std::string> revfiles;
	split(forward_fastqs,delim,forfiles);
	split(rev_fastqs,delim,revfiles);
	
	if (revfiles[0].compare("null")==0 | forfiles[0].compare("null")==0) {
	  std::cerr << "Missing forward or reverse fastq file and program is in paired-end mode\n";
	  exit(1);
	}

      	if (revfiles.size() > forfiles.size()) {std::cerr << "Incorrect number of forward/reverse read files, too few forward read files provided\n";}
	for (int ii=0;ii<forfiles.size();++ii) {
	  opt.files.push_back(forfiles[ii]);
	  opt.file_direction.push_back(true);
	  std::cerr << revfiles[ii] << " SIZE\n";
	  if (ii < revfiles.size()) {
	    opt.files.push_back(revfiles[ii]);
	    opt.file_direction.push_back(false);
	  } else {
	    std::cerr << "Incorrect number of forward/reverse read files, too few reverse read files provided\n";
	  }
	}
      } else {
	split(forward_fastqs,delim,opt.files);
      }

      std::cout.sync_with_stdio(false);
      setvbuf(stdout, NULL, _IOFBF, 1048576);

      if (opt.files[0].compare("null")==0) {
	if (command=="quantify") {
	  std::cerr << "Error, missing fastq file. Enter using -f\n";
	  exit(1);
	}
      } else {
	for (int ii=0;ii<opt.files.size();++ii) {	    
	  std::ifstream test;
	  test.open(opt.files[ii].c_str(), std::ios::in);
	  if (!test.is_open()) {
	    std::cerr << "Unable to open fastq file " << opt.files[ii] << std::endl;
	    exit(1);
	  }
	  test.close();
	}
      }

      if (opt.taxonomy_files[0].compare("null")==0) {
	if (command=="quantify") {
	  std::cerr << "Error, missing taxonomy file. Enter using -t\n";
	  exit(1);
	}
      } else {
	for (int ii=0;ii<opt.taxonomy_files.size();++ii) {	    
	  std::ifstream test;
	  test.open(opt.taxonomy_files[ii].c_str(), std::ios::in);
	  if (!test.is_open()) {
	    std::cerr << "Unable to open taxonomy file " << opt.taxonomy_files[ii] << std::endl;
	    exit(1);
	  }
	  test.close();
	}
      }

    if (command=="index") {

      if(checkIndexOptions(opt)) {

	Kmer::set_k(opt.k);
	KmerIndex index(opt);
	index.BuildTranscripts(opt);
	index.write(opt.index);

      } else {
	exit(1);
      }

    }

    if (command=="quantify") {

      
      std::stringstream logfile_track;

      std::cerr << "Reading fasta index\n";
      logfile_track << "Reading fasta index\n";
      fastaIndex findex(opt.transfasta);

      LKFilter like_filter(opt);
      if (opt.harp_filter) {
	std::cerr << "Calculating likelihood based read filter parameters\n";
	logfile_track <<  "Calculating likelihood based read filter parameters\n";
	like_filter.GetParameterValues(opt);
	like_filter.use = true;
	std::cerr << "Likelihood filter parameters: mean=" << like_filter.mean << " SD=" << like_filter.sd << std::endl;
	logfile_track << "Likelihood filter parameters: mean=" << like_filter.mean << " SD=" << like_filter.sd << std::endl;
     
      } else {
	std::cerr << "Not performing likelihood based filter of reads\n";
	logfile_track << "Not performing likelihood based filter of reads\n";	
      }
      
      HLK likelihoods(findex,opt,logfile_track);
      double number_mapped;
 
      EarlyTaxonomy earlytax;
      earlytax.BuildEarlyTaxonomy(opt.taxonomy_files,findex);


	std::chrono::system_clock::time_point tp;
	tp = std::chrono::system_clock::now();
      
	{
	  KmerIndex index(opt);
	  index.load(opt,logfile_track);
	  
	  MinCollector collection(index, opt);
	  unsigned int num_processed = ProcessReads(index,opt,collection,findex,logfile_track,like_filter,earlytax,likelihoods);
	  likelihoods.nummapped = (double) num_processed;

	  if (opt.minimum_frequency_cutoff < 0) {
	    std::cout << "Setting frequency to 1/" << likelihoods.nummapped << std::endl;
	    opt.minimum_frequency_cutoff = 1.0/likelihoods.nummapped;
	    likelihoods.minimum_frequency_cutoff = opt.minimum_frequency_cutoff;
	  }

	}
	
	auto diff = std::chrono::system_clock::now() - tp;
	tp = std::chrono::system_clock::now();
	std::cout << "Pseudomapping/aligning:  " << std::chrono::duration_cast<chrono::seconds>(diff).count() << " sec" << std::endl;
 
	likelihoods.estimate_haplotype_frequencies();

	diff = std::chrono::system_clock::now() - tp;
	tp = std::chrono::system_clock::now();     
	std::cout << "EM:  " << std::chrono::duration_cast<chrono::seconds>(diff).count() << " sec" << std::endl;

	std::string outer = opt.out + ".freqs";
	std::ofstream outfile;
	outfile.open(outer.c_str(), std::ios::out);
	if (!outfile.is_open()) {std::cerr << "Unable to open output file\n";exit(1);}
	
	if (opt.collapse==true) {
	  earlytax.output(likelihoods,outfile);
	} else {
	  earlytax.output2(likelihoods,outfile);
	}

	std::string log = opt.out + ".log";
	std::ofstream logfile;
	logfile.open(log.c_str(), std::ios::out);
	if (!logfile.is_open()) {std::cerr << "Unable to open log file\n";exit(1);}
	printConfig(logfile,opt);
	logfile << logfile_track.str();
	logfile.close();
    }

    if (command=="tabulate") {

      if (all_samples.compare("null")==0) {
	std::cerr << "Command \"tabulate\" without input sample filelist. Enter with --samples\n";
	exit(1);
      } else {
	ResTable table(opt);
	table.tabulate(all_samples);
	table.output();
 
      }
    }   
    

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
  return 0;
}

    
    
