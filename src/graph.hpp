#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <math.h>

#define MAXN 1000
// max graph size to avoid numerical errors

class Graph{
public:
    int n;
    std::vector<std::string> node_attr;
    std::vector<std::vector<std::pair<int, std::string>>> adj;
    int ni = 0;
};



#endif
