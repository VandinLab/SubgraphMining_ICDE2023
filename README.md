# VC-dimension and Rademacher Averages of Subgraphs, with Applications to Graph Mining

This is the code and data repository for the paper **_VC-dimension and Rademacher Averages of Subgraphs, with Applications to Graph Mining_**, presented at [IEEE ICDE 2023](https://icde2023.ics.uci.edu).


### Compilation
The code depends on the [NLopt](https://nlopt.readthedocs.io/en/latest/) and the [VF3Lib](https://github.com/MiviaLab/vf3lib) libraries, 
which are provided in this repository for convenience. 
The file ```init.sh``` provides a script to compile and install the dependencies locally, and to complile the source code for our tools. 
The executable of our tools will be located in the ```src```  subdirectory.

The code is tested only with `gcc-11` under Linux and macOS (Intel) operating systems, mut might run also with other C++ compilers and operating systems.

To execute the scrupt, simply run:

    ./init.sh

### Frequent subgraph mining
The first tool, ```frequent_subgraph_sampler```, given a dataset of graphs and parameters $\epsilon, \delta > 0$, returns a sample of the original
dataset such that running a frequent subgraph miner on the sample yields, with probability $1-\delta$, an $\epsilon$-approximation of the 
frequent subgraphs that would be obtained by running the miner on the entire dataset. 

The program takes the following inputs:

    ./frequent_subgraph_sampler [-m size] [-s seed] [-r] [-l] eps delta n input output
    
with:
- `-m size`, the maximum number of nodes in the subgraphs of interest (optional, default = $\infty$)
- `-s seed`, the seed (optional, default = 0)
- `-r`, if set, the tool uses Rademacher averages, otherwise it uses the VC dimension 
- `-l`, if set, the tool takes into account node labels
- `eps`, the approximation quality
- `delta`, probability of having any false negatives
- `n`, the number of graphs in the original dataset
- `input`, path to the input dataset file
- `output`, path for the output dataset file

### True frequent subgraph mining
The second tool, ```true_frequent_subgraph_evaluator```, given a dataset of graphs and a parameter $\delta > 0$, returns a value $\epsilon$ such that, 
with probability $1-\delta$, the maximum deviation of the frequency of any subgraph in the dataset from its true frequency in the underlying generative 
distribution is bounded by $\epsilon$.

The program takes the following inputs:

    ./true_frequent_subgraph_evaluator [-m size] [-s seed] [-r] [-l] delta n input 
    
with:
- `-m size`, the maximum number of nodes in the subgraphs of interest (optional, default = $\infty$)
- `-s seed`, the seed (optional, default = 0)
- `-r`, if set, the tool uses Rademacher averages, otherwise it uses the VC dimension 
- `-l`, if set, the tool takes into account node labels
- `delta`, probability of error
- `n`, the number of graphs in the original dataset
- `input`, path to the input dataset file

### File formats
Our tools use the graph format used by most subgraph miners:

    t # <graph-id>
    v <vertex-id> <vertex-label>
    ...
    e <vertex-id> <vertex-id> <edge-label>
    ...
    
with graph and vertex ids contiguous integers starting from 0, and vertex and node labels interers.


