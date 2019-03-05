# Haystack

This software implements the analytical cache model discussed by the paper "A Fast Analytical Model of Fully Associative Caches". The software was developed by SPCL (ETH Zurich).

## Installation

Before installing the package make sure you install die following dependencies:
- clang
- gmp
- ntl
- boost (program options)
- libyaml

On Ubuntu systems (18.04) you can install these dependencies using the following command:
```
sudo apt-get install libclang-dev libgmp3-dev libntl-dev libboost-program-options-dev libyaml-dev
```
Once all dependencies are installed change to the haystack folder and run the commands:
```
./get_submodules.sh
./autogen.sh
```
The two commands get all submodules and initialize autotools. After these steps we are ready to configure and build the project:
```
./configure --prefix=$HOME
make
make install
```
We set the prefix to the home folder to install the tool in the user directory (this is optional). Additional configure options also allow us to point autotools to dependencies such as the boost program options library (in case they are not found automatically).

## Usage

Haystack is a command line tool that analyzes c source files with an annotated scop. The example folder contains two valid test input files:
- gemm.c
- gemm.tiled.c

All other files are used by the unit tests and due to undefined loop bounds don't work out of the box. To analyze the cache misses of gemm, we type the following command:
```
haystack -f ../examples/gemm.c
```
The tool then reports the number of cache misses per statement/memory reference:
```
...
12     for (int k = 0; k < D; k++) 
13       for (int j = 0; j < D; j++)
14         C[i][j] += alpha * A[i][k] * B[k][j];
   -------------------------------------------------------------------
             ref  type  comp[%]  L1[%]    L2[%]    tot[%]   reuse[ln]
         C[i][j]  rd    0.0000   0.0000   0.0000   24.9878  11,14
         A[i][k]  rd    0.0015   0.0000   0.0000   24.9878  14
         B[k][j]  rd    0.0015   1.5602   1.5602   24.9878  14
         C[i][j]  wr    0.0000   0.0000   0.0000   24.9878  14
...
```
The columns provide the following information:
- ref: memory reference
- type: read or write access
- comp: compulsory misses in percent of the total number of memory accesses
- L1: capacity misses (L1) in percent of the total number of memory accesses
- L2: capacity misses (L2) in percent of the total number of memory accesses
- tot: number of memory accesses in percent of the total number of memory accesses
- reuse: line numbers that contain the memory references that last accessed the same cache line (reuse)

The tool additionally reports the absolute numbers of cache misses for the entire scop:
```
compulsory:                  196'608
capacity (L1)             67'043'328
capacity (L2)             67'043'328
total:                 4'297'064'448
```
The tool provides a number of additional program options:
```
haystack -h
Program options:
  -h [ --help ]                         print the program options
  -c [ --cache-sizes ] arg (=32768 524288)
                                        cache sizes in byte
  -l [ --line-size ] arg (=64)          cache-line size in byte
  -f [ --input-file ] arg               set the source file [file name]
  -I [ --include-path ] arg             set the include path [include path]
  -s [ --scop-function ] arg            set the scop function scop
```
We can use the -c and -l options to define different cache configurations and the -I option to specify one or more include paths. The -s option allows us to select the function from which we want to extract the scop. This is useful for input files with multiple scops.
