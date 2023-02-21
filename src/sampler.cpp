#include <vector>
#include <string>
#include "subgraph_sampling.hpp"
#include "utils.hpp"

using namespace std;


int main(int argc, char *argv[]){
    if ( argc < 6 || argc > 6 ) {
      cerr << "Parameters: n input to_sample output seed" << endl;
      return 1;
    }

    int n = atoi(argv[1]);
    string in = argv[2];
    int to_sample = atoi(argv[3]);
    string out = argv[4];
    int seed = atoi(argv[5]);

    ifstream file(in);
    ofstream ofile;
    ofile.open(out);

    std::mt19937 mt = std::mt19937( seed >= 0 ? seed : time(0));

    vector<int> sampled_idxs = vector<int>();
    while(sampled_idxs.size() < to_sample){
        int s = mt()%n;
        sampled_idxs.push_back(s);
    }
    sort(sampled_idxs.begin(), sampled_idxs.end());


    FILE* cfile = fopen(in.c_str(), "r");
    int id = 0, cnt = 0;
    char array[100];
    char* err = fgets ( array, 100, cfile );
    while ( !feof ( cfile ) && id < sampled_idxs.size() ) {
        Graph tau; readGraph(cfile, tau, true);

        while(cnt == sampled_idxs[id]){
            ofile << "t # " << cnt << endl;
            for(int i=0; i<tau.n; i++)
                ofile << "v " << i << " " << tau.node_attr[i] << endl;
            for(int i=0; i<tau.n; i++){
                for(int j=0; j<tau.adj[i].size(); j++){
                    if(tau.adj[i][j].first > i)
                    ofile << "e " << i << " "<< tau.adj[i][j].first << " " << tau.adj[i][j].second << endl;
                }
            }
            id++;
        }
        cnt++;
    }
    cerr << "Written " << id << " transactions" << endl;
    ofile.close();
}

