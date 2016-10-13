# Karp

Karp is a program for identifying the relative frequencies of reference sequences contributing to a pooled DNA sample. Karp was developed with 16S rRNA sequencing as its primary application, but will work in any context where reference sequences are available in fasta format and the pooled DNA is in fastq format with base quality scores. Karp employs a pseudoalignment step to quickly match reads with potential references. Reads are locally aligned to the references they pseudoalign with, and Karp calculates a likelihood of each read originating from each reference using base quality scores. Finally, Karp employs an EM algorithm that uses information from all the reads to accurately estimate the relative frequencies of each reference in the sample.

## Installation

After cloning the repository, or downloading and decompressing the tarball follow the following steps to install Karp. In the Karp directory
      
      mkdir build
      cd build
      cmake ..
      make

This should create a working karp executable in the folder /src

###Trouble Shooting

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

 


