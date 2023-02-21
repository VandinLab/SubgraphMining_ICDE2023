#include <vector>
#include <string>
#include "subgraph_sampling.hpp"

using namespace std;

int nextTransaction(std::ifstream& file, Graph& tau, bool labeled){
    using namespace std;
    bool read = true;
    string line, tmp, tmp2, tmp3;
    getline(file, line);
    getline(file, line);
    while(line[0] == 'v'){
        if(line.back() == '\r') line.pop_back();
        auto ss = stringstream(line);
        getline(ss, tmp, ' '); getline(ss, tmp, ' '); getline(ss, tmp, ' ');
        if(!labeled) tmp = "0";
        tau.node_attr.push_back(tmp);
        getline(file, line);
    }
    tau.n = tau.node_attr.size();
    tau.adj = vector<vector<pair<int, string>>>(tau.n, vector<pair<int, string>>(0));
    while(read && line[0] == 'e'){
        if(line.back() == '\r') line.pop_back();
        auto ss = stringstream(line);
        getline(ss, tmp, ' ');getline(ss, tmp, ' ');getline(ss, tmp2, ' ');getline(ss, tmp3, ' ');
        if(!labeled) tmp3 = "0";
        tau.adj[stoi(tmp)].push_back(make_pair(stoi(tmp2), tmp3));
        tau.adj[stoi(tmp2)].push_back(make_pair(stoi(tmp), tmp3));
        read = (bool)getline(file, line);
    }

    if(!read || line[0]=='\n') return -1;
    return stoi(line.substr(2, line.size()-2));
}



int main(int argc, char *argv[]){
    if(argc - optind != 4){
        cerr << "Parameters: complete_input n samp_input n  Provided " << argc << " parameters" << endl;
        return 1;
    }
    string in1 = argv[1];
    string in2 = argv[3];

    ifstream file(in1);
    ifstream file2(in2);
    int nc = atoi(argv[2]), ns = atoi(argv[4]);

    vector<Graph> complp, samplp;
    int cntc = 0, cnts = 0;
    string line;
    getline(file, line);
    int fr = stoi(line.substr(2, line.size()-2));
    while (fr > -1) {
        Graph tau; tau.ni = fr;
        fr = nextTransaction(file, tau, 1);
        cntc++;
        complp.push_back(tau);
    }
    file.close();

    getline(file2, line);
    fr = stoi(line.substr(2, line.size()-2));
    while (fr > -1) {
        Graph tau; tau.ni = fr;
        fr = nextTransaction(file2, tau, 1);
        cnts++;
        samplp.push_back(tau);
    }
    file2.close();

    int cnt_matches = 0;
    for(int i=0; i<samplp.size(); i++){
        bool ok = 0;
        for(int j=0; j<complp.size(); j++){
            Matcher m = Matcher(samplp[i], complp[j]);
            if(m.isIsomorphic()){
                cout << i << " " << j << " "<<((double)samplp[i].ni / ns) << " "<< ((double)complp[j].ni / nc)<<" " << ((double)samplp[i].ni / ns) - ((double)complp[j].ni / nc) << endl;
                ok = 1;
                cnt_matches++;
                break;
            }
        }
        if(!ok){
            //cout << i << " with freq. " << ((double)samplp[i].ni / ns) << "is not freq. in the complete dataset" << endl;
            cout << i << " -1 " << ((double)samplp[i].ni / ns) << " 0.0 0.0" << endl;
        }
    }

    if(complp.size() - cnt_matches > 0)
        cout << complp.size() - cnt_matches <<  "freq. patterns not freq. in the sampled dataset" << endl;

}
