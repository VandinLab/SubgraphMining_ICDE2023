#ifndef VC_H
#define VC_H

#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <set>
#include "graph.hpp"
#include "vf3_wrappers.hpp"
#include "utils.hpp"

#define VC_CONST 0.5
#define CHECK_ISO 0

struct VCSampler{
public:
    VCSampler(std::string _file, bool _lab, int seed = -1){
        using namespace std;
        infile = _file;
        labeled = _lab;
        mt = std::mt19937( (seed >= 0 ? seed : time(0) ) );

        C = vector<vector<ld>>(MAXN, vector<ld>(MAXN));
        // precompute binomial coeffcient
        for (int i = 0; i < MAXN; i++) {
            for (int j = 0; j < min(i, MAXN); j++) {
                if (j == 0 || j == i)
                    C[i][j] = 1.0;
                else
                    C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }

    }

    int compute_vc_and_sample(double eps, double delta, int k = -1){
        using namespace std;
        int vc;
        if(k < 0 || k > MAX_P_SIZE)
            vc = get_d_bound();
        else
            vc = get_c_bound(k);

        int to_sample = min( n, (int)ceil(VC_CONST/(eps*eps)*(vc + log(1/delta))) );
        if(VERBOSE) cerr << "VC: " << vc << " samples:" << VC_CONST/(eps*eps)*(vc + log(1/delta)) << " n: " << n <<endl;
        sampled_idxs = vector<int>();
        while(sampled_idxs.size() < to_sample){
            int s = mt()%n;
            sampled_idxs.push_back(s);
        }
        return vc;
    }


    int output_sampled_dataset(std::string filename){
        using namespace std;
        FILE* cfile = fopen(infile.c_str(), "r");
        ofstream ofile;
        ofile.open(filename);
        sort(sampled_idxs.begin(), sampled_idxs.end());

        int cnt = 0, id = 0;
        char array[100];
        char* err = fgets ( array, 100, cfile );
        while ( !feof ( cfile ) && id < sampled_idxs.size() ) {
            Graph tau;
            if(cnt == sampled_idxs[id]) readGraph(cfile, tau, true);
            else skipGraph(cfile); // faster

            while(cnt == sampled_idxs[id]){
                ofile << "t # " << cnt << "\n";
                for(int i=0; i<tau.n; i++)
                    ofile << "v " << i << " " << tau.node_attr[i] << "\n";
                for(int i=0; i<tau.n; i++){
                    for(int j=0; j<tau.adj[i].size(); j++){
                        if(tau.adj[i][j].first > i)
                        ofile << "e " << i << " "<< tau.adj[i][j].first << " " << tau.adj[i][j].second << "\n";
                    }
                }

                id++;
            }

            cnt++;
        }
        fclose(cfile);
        ofile.close();
        if(VERBOSE) cerr << "Written " << id << " transactions" << endl;
        return id;
    }


    int get_d_bound(){
        using namespace std;
        FILE* cfile = fopen(infile.c_str(), "r");
        vector<Graph> T;
        int q = 0, cnt = 0;
        char array[100];
        char* err = fgets ( array, 100, cfile );
        while ( !feof ( cfile ) ) {
            Graph tau; readGraph(cfile, tau, true);
            cnt++;
            if(tau.n <= q) continue;
            bool ok = true;
            for(int i=0; CHECK_ISO && i<T.size(); i++){
                if(T[i].n != tau.n) continue;
                Matcher m = Matcher(T[i], tau);
                if(m.isIsomorphic()) {ok = false; break;}
            }
            if(!ok) continue;

            T.push_back(tau);
            int minsize = 1000000000, mini = -1;
            for(int i=0; i<T.size(); i++){
                if(T[i].n < minsize){
                    minsize = T[i].n;
                    mini = i;
                }
            }
            if(minsize > q){
                q++;
            }
            else{
                T[mini] = T[q];
                T.pop_back();
            }
        }
        n = cnt;
        fclose(cfile);
        return q;
    }

    int get_c_bound(int k){
        using namespace std;
        FILE* cfile = fopen(infile.c_str(), "r");
        vector<Graph> T;
        int q = 0, cnt = 0;
        char array[100];
        char* err = fgets ( array, 100, cfile );

        while ( !feof ( cfile ) ) {
            Graph tau; readGraph(cfile, tau, true);
            cnt++;
            if(tau.n <= q) continue;
            long double ni = log2((long double)tau.n); // n choose 1
            for(int j=2; j <= min(k, tau.n); j++){
                long double log_binom = 0.0;
                if(tau.n < MAXN){
                    log_binom = log2(C[tau.n][j]);
                } else {
                    for(int j2 = 0; j2 < j; j2++) log_binom += log2((long double)(tau.n-j2));
                    for(int j2 = 1; j2 <= j; j2++) log_binom -= log2((long double)j2);
                }
                if(!labeled && j < 15 && log2(n_connected_graphs[j]) < log_binom)
                    ni = logsumexp2(ni, log2(n_connected_graphs[j]));
                else
                    ni = logsumexp2(ni, log_binom);
            }
            tau.ni = (int)floor(ni) + 1;
            if(tau.ni <= q) continue;
            bool ok = true;
            for(int i=0; CHECK_ISO && i<T.size(); i++){
                if(T[i].n != tau.n) continue;
                Matcher m = Matcher(T[i], tau);
                if(m.isIsomorphic()) {ok = false; break;} // could test also subgraph isomorphism to help avoid chains, bit it is slower
            }
            if(!ok) continue;

            T.push_back(tau);
            int minsize = 1000000000, mini = -1;
            for(int i=0; i<T.size(); i++){
                if(T[i].ni < minsize){
                    minsize = T[i].ni;
                    mini = i;
                }
            }
            if(minsize > q){
                q++;
            }
            else{
                T[mini] = T[q];
                T.pop_back();
            }
        }
        n = cnt;
        fclose(cfile);
        return q;
    }


    double get_epsilon_evc(double delta, int k){
        int evc;
        if(k < 0 || k > MAX_P_SIZE)
            evc = get_d_bound();
        else
            evc = get_c_bound(k);
        if(VERBOSE) std::cerr << "Evc: "<<evc << " n: " << n<< std::endl;
        double eta = 2.0*sqrt(2.0*evc*log(n+1)/n) + sqrt(2.0*log(2.0/delta)/n);
        return eta;
    }


    int n, labeled;
    std::mt19937 mt;
    std::string infile;
    std::vector<int> sampled_idxs;
    std::vector<std::vector<ld>> C;
    const double n_connected_graphs[16] = {1.0, 1.0, 1.0, 2.0, 6.0, 21.0, 112.0, 853.0, 11117.0, 261080.0, 11716571.0, 1006700565.0, 164059830476.0, 50335907869219.0, 29003487462848061.0, 31397381142761241960.0};
};



#endif
