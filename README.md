# Karp

Karp is a program for identifying the relative frequencies of reference sequences contributing to a pooled DNA sample. Karp was developed with 16S rRNA sequencing as its primary application, but will work in any context where reference sequences are available in fasta format and the pooled DNA is in fastq format with base quality scores. Karp employs a pseudoalignment step to quickly match reads with potential references. Reads are locally aligned to the references they pseudoalign with, and Karp calculates a likelihood of each read originating from each reference using base quality scores. Finally, Karp employs an EM algorithm that uses information from all the reads to accurately estimate the relative frequencies of each reference in the sample.

## Installation

After cloning the repository, or downloading and decompressing the tarball follow the following steps to install Karp. In the Karp directory
      
      mkdir build
      cd build
      cmake ..
      make

This should create a working karp executable in the folder /src
If you get an error message during compilation see the bottom of this readme for trouble shooting tips.

## Usage

Command    |    Description
------------------------ | -------------------------------------------------------------------------------
**-c / -\-command**    | Which of Karp's functions should run, options are 'index' 'quantify' 'tabulate'
**-h / -\-help**       | Print out Karp options and their descriptions

#### Indexing
The first stage of analysis with Karp is building the k-mer index for pseudoalignment. Indexing is performed if the option 'index' is specified as the **-c** value.

Command    |         Description
------------------ | -------------------------------------------------------------------------------------
**-r / -\-ref**        | Reference fasta file to build index from, must have matching .idx file created with samtools faidx.<br> Enter multiple files with comma delimiter.
**-i / -\-index**  | Name of index file to output
**-k / -\-kmer**   | Length of k-mer to use when building index, must be odd and range between 3 and 31 [default = 31].


#### Quantifying
The main function of Karp is to quantify the taxonomy in a pooled DNA sample, that is done after the index has been constructed using the option 'quantify' as the **-c** value.

 Command                              |    Description
---------------------------------------- | --------------------------------
**-r / -\-ref**          | Reference fasta files, must have matching .idx file created with samtools faidx.<br>Multiple files entered with comma delimiter.
**-i / -\-index**        | Name of k-mer index file built with 'index'
**-f / -\-forward**      | Fastq files to be quantified, can be gzipped. Enter multiple files with comma<br>separating them. If quantifying single-end reads, enter fastq files with<br>this command. If quantifying paired-end reads enter forward files with <br>this command.
**-\-paired**            | Enter this flag if you are quantifying paired-end reads
**-q / -\-reverse**      | Reverse oriented fastq files to be quantified. Should match order of files<br>entered with -\-forward. Can be gzipped.
**-t / -\-tax** | Taxonomy file(s) matching reference fasta file(s). See note below for more<br>details about acceptable formats.
**-o / -\-out** | Base name for output files. Karp will output a '.freqs' file containing results<br>and '.log' file [default = 'karp'].
**-\-threads**        | Number of threads to use [default = 1].
**-\-phred**          | Version of phred quality scores being used. Default is Phred+33,<br>used by Illumina 1.8+. Other valid option is '64' which corresponds<br>to Phred+64 quality coding [default = 33].
**-\-min\_freq**      | Lower frequency bound during EM update step. Taxa with frequencies below<br>this threshold are removed from the analysis [default = 1/Number of reads].
**-\-like\_thresh**   | The maximum likelihood z-score threshold. Karp uses distribution of quality<br>scores in input fastq file to calculate a "null" distribution for likelihood<br>values when each read is matched to the correct reference and all mismatches<br>are the result of sequencing error. Using this distribution Karp calculates<br>a z-score for the maximum likelihood observed for each read, and if this<br>value falls to far outside the theoretical distribution the read is filtered<br>from analysis [default = -2.0].
**-\-collapse**       | Run Karp in Collapse mode. In collapse mode classification is done at<br>the level of taxonomic labels rather than individual reference sequences.


##### Advanced Options for Quantifying

Command    |    Description
 --------------------- |   ------------
**-\-max\_like\_out**    | Output a file containing the maximum likelihood value for each read,<br>use this to diagnose appropriate max likelihood filter threshold.
**-\-strict**            | Use the strict intersection of equivalence classes to declare psuedoalignment of read<br>to reference (this is the definition of pseudoalignment used by Kallisto).<br>Default Karp behavior is that if no strict intersections exists, the reference(s) with<br>the greatest number of matching k-mers are declared as pseudoaligning to a query read. 
**-\-em\_tolerance**    | When the squared sum of allele frequency changes between EM algorithm iterations<br>falls below this threshold the algorithm has converged [default = 1e-12].
**-\-no\_harp\_filter**  | Turn off the likelihood z-score filter.
**-\-fail**            | Output a file listing the reads that fail to pseudoalign and those that are excluded by<br>the maximum likelihood z-score filter.
**-\-max\_em\_it**       | Maximum number of iterations of EM algorithm before declaring failure<br>to converge [default = 1000].

#### Tabulate
After multiple samples have been classified using Karp, the program can be used to calculate some basic summary statistics using the 'tabulate' command.

Command    |    Description
--------------------- |   ------------
**-\-samples**           | File with list of Karp output files from each sample that has been classified. File should have two columns,<br>the first with a sample ID and the second with a path to the corresponding Karp output file

## Best Practices

If using Karp with multiple threads, the internal local aligner runs more quickly if the reference database fasta file is broken into multiple smaller files. A single kmer index can still be used, as long as the references in the multiple fasta files are all present in the single index.

The default likelihood threshold that karp uses to filter reads for being unlikely to originate from the reference database can be too lenient or strict with some distributions of read quality scores. If it seems that karp is filtering too many (or too few) reads, it can be helpful to use the ** -\-max\_likes\_out** option, and plot the maximum likelihood for each read to get a sense of how to adjust the ** -\-like\_thresh** option.

Data from a single 16S hypervariable region often have little power to distinguish between closely related reference sequences. When this is the case Karp will partition the probability weight across the closely related species, and they will appear at low frequency. In such situations the default minimum EM update frequency can be too strict, and setting it to a lower value based on the depth of the samples being classified can improve results. The default value is 1/Number of Reads, with a single hypervariable region 1/(10\*Number of Reads) or 1/(20\*Number of Reads) can sometimes give more accurate results.

## Test Example

In the repo there is a directory labeled `example`, this contains some toy files for demonstrating how to use Karp.

| File                | Description                                                |
| -------------------------------- | ---------------------------------------------------------- |
| reference.fasta                  | Reference database in fasta format, contains 500 sequences |
| reference.fasta.fai              | Tabix index of reference fasta file                        |
| reference.tax                    | Taxonomy file with description of sequences in reference.fasta   |
| simulated.data.fastq.gz          | Fastq file containing 10,000 reads with mixture of organisims from reference.fasta |
| simulated.data.actual_counts.txt | True origins of reads in simulated.data.fastq.gz |

The first step of running Karp is to create a k-mer index for pseudoalignment:

``
./karp -c index -r reference.fasta -i reference.index
``  
  
This will create the index file `reference.index`. Next we quantify the simulated reads:

``
./karp -c quantify -r reference.fasta -i reference.index -f simulated.data.fastq.gz -o simulated.results -t reference.tax 
``

In the example folder we now have the files `simulated.results.freqs` and `simulated.results.log`. `simulated.results.freqs` contains Karp's estimates of the number of reads each reference sequence in `reference.fasta` contributed to the simulated sample. `simulated.results.log` contatins information about the Karp run we just performed.

## Trouble Shooting

Error message: `unrecognized command line option "-std=c++0x"`

or 

Error message: `error: call of overloaded ‘to_string(size_t&)’ is ambiguous`

Likely your default C++ compiler is out of date, check for newer compiler availability. Once you have found/installed an up to date compiler, delete and recreate the build directory and rerun the cmake step with command

`cmake -DCMAKE_CXX_COMPILER=/path/to/compiler -DCMAKE_C_COMPILER=/path/to/compiler ..`

 Error message: `Could NOT find HDF5 (missing: HDF5_LIBRARIES HDF5_INCLUDE_DIRS
  HDF5_HL_LIBRARIES)`

 Check if you have HDF5 libraries installed on your machine. If you do not they can be downloaded and installed [https://support.hdfgroup.org/HDF5/](https://support.hdfgroup.org/HDF5/ "HDF5 Homepage")

 If you have HDF5 installed, and are still getting this message, make sure your $HDF5_DIR variable is pointing to the correct directory. If it is, try adding -DCMAKE_PREFIX_PATH=/path/to/HDF5 to the cmake command.

 There is also a known issue with some versions of CMAKE and the HDF5 library. If you are running CMAKE versions 3.6.0 or 3.6.1 see if you have access to an earlier or later version of CMAKE and try using it instead.

 


