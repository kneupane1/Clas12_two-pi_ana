# Clas12_sim_checks

This code is designed for the two-pion channel analysis using CLAS12 RGA data, which has been converted into ROOT files through dst2root. This program is build on c++/root. It is a modified version of Nick Tyler's CLAS12 analysis framework, tailored specifically for the two-pion channel and extended to include simulation analysis.

Prerequisites:

- [cern root](https://root.cern.ch/)

### cpp

To build:

1. clone this code or download from github : git clone https://github.com/kneupane1/Clas12_two-pi_ana.git
2. cd Clas12_two-pi_ana
3. mkdir build
4. cd build
5. cmake ..
6. make
7. To run:

   ./clas12_analysis output.root /path/to/input/input_file.root

   OR to use multi thrades:
   CLAS12_E=10.6041 NUM_THREADS=4 ./clas12_analysis output.root /path/to/input/input_file.root

To avoid issues, ensure each thread processes at least 2 files, and the number of threads does not exceed the number of cores. For example, with 16 files on a 4-core system, set NUM_THREADS=4 (4 files per thread). For 4 files on the same system, use NUM_THREADS=1 or 2 (4 or 2 files per thread, respectively).
