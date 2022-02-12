# CoRBeC++: An Exact Solver for Cluster Editing Problem

The **Co**nflict **R**educed **Be**st-fit **C**luster is a software for
solving the cluster editing problem. The cluster editing problem is to transform an input graph into a cluster graph (a disjoint union of complete graphs) by performing a minimum number of edge-modifications (addition or deletion). 
This repository provides an exact solver for the cluster editing problem.
For a detailed overview of the techniques used in our solver, please refer to the following resource:

- Exact Solver ([PDF](exact_description))

![cluster-editing](https://user-images.githubusercontent.com/9654047/119774492-88069e00-bec2-11eb-8800-c4abfcacb82f.png)

Requirements
-----------

 - A 64-bit Linux operating system
 - A modern, ![C++17](https://img.shields.io/badge/C++-17-blue.svg?style=flat) ready compiler
 - Gurobi ILP Solver (Version >= 9.5) ([download](https://www.gurobi.com/downloads/))

Build Application
-----------

1. Clone the repository including submodules:

   ```git clone --recursive https://github.com/sachin-4099/Cluster-Editing-Exact-Solver```
2. Build the binary:

   If the Gurobi Version is 9.5 (make changes in the command depending on the version number):
   
   ```g++ -O2 corbec++.cpp -o corbec++ -I/opt/gurobi950/linux64/include -L/opt/gurobi950/linux64/lib -lgurobi_c++ -lgurobi95```

Run Application
-----------

Reads a cluster editing instance (a graph file in `.gr` format) from stdin and prints the modifications along with the details of steps performed to stdout.

A usage example would be:

    cat tests/exact011.gr | ./corbec++ > modifications.txt

By default, the solver runs for `30 minutes` and the number of threads is set to `4`.

The time limit of the solver can be adjusted via the `--time-limit` flag (in seconds):

    cat tests/exact011.gr | ./corbec++ --time-limit=100 > modifications.txt

The number of threads used by the solver can be adjusted via the `--num-threads` flag:

    cat tests/exact011.gr | ./corbec++ --num-threads=56 > modifications.txt
