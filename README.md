# Promatch

Promatch is an adaptive locality-aware greedy predecoder that extends real-time quantum error correction for surface codes of distances 11 and 13 ([click here for the paper](https://dl.acm.org/doi/abs/10.1145/3620666.3651339)). This repository includes source code for simulating and evaluating our proposed designs in Promatch. It focuses on very low logical error rates using rare event simulation. You can use this repository to reproduce important findings from the paper, specifically in Section 6.1 (Table 2 and Table 3 on Logical Error Rates) and Section 6.2 (Figure 14 and Figure 15 on Error Rate Sensitivity Analysis up to 5×10−4).

## Hardware dependencies
Any computing clusters should be capable of executing the necessary experiments. For our evaluations, we utilized 120 cores. It is recommended to allocate 2GB of memory for each core.

## Software dependencies:
The code uses CMake version 3.20.2, but older versions also work with adjustments to the CMakeLists.txt file. It is compiled with the g++ compiler supporting C++17 and requires OpenMPI version 4.x.x. All additional required software is included with the code, and CMake manages them. This code has been tested on MacOS and Linux.

## Installation:
The commands below create the promatch executable to run experiments.
```
cd Release
cmake .. -DCMAKE\_BUILD\_TYPE=Release
make
```
