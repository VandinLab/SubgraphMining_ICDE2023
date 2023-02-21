#ifndef RADEMACHER_H
#define RADEMACHER_H
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <math.h>
#include <random>
#include <algorithm>
#include "graph.hpp"
#include "vf3_wrappers.hpp"
#include "utils.hpp"
#include "mynlopt.h"

#define ll long long
#define ld long double



static double opt_fun(unsigned nn, const double* x, double* grad, void* datap){
    using namespace std;
    vector<pair<long double, long double>> data = *((vector<pair<long double, long double>>*) datap);
    vector<long double> exponents(data.size());
    long double s_square = 1.0L*x[0]*x[0];
    for(int i=0; i<data.size(); i++){
        exponents[i] = data[i].first + data[i].second * s_square;
    }
    long double log_sum = logsumexp(exponents);
    return (double)(log_sum/x[0]);
}



struct LabeledRademacherComputer{
public:
    LabeledRademacherComputer(std::string _file, int _n, int _maxs, int seed, int _maxdepth){
        using namespace std;
        infile = _file;
        n = _n;
        maxdepth = _maxdepth;
        mt = std::mt19937( (seed >= 0 ? seed : time(0) ) );
        if(_maxs < 0 || _maxs > MAX_P_SIZE) maxs = MAX_P_SIZE;
        else maxs = _maxs;

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

    double get_epsilon(double delta){
        double r_tilde = sample_and_update(n);
        double eta = 2.0*r_tilde + sqrt(2.0*log(2.0/delta)/n);
        return eta;
    }



    void progressive_sampling(double eps, double delta){
        int samp_size = std::min(n, (int)(8.0*log(2.0/delta)/(eps*eps)) );
        while(true){
            double r_tilde = sample_and_update(samp_size);
            double eta = 2.0*r_tilde + sqrt(2.0*log(2.0/delta)/samp_size);
            if(VERBOSE) std::cerr << "Eta " << eta << std::endl;
            if(eta < eps || samp_size == n){
                break;
            }
            samp_size = next_sample_size(eta, eps, samp_size);
        }
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



private:

    void count_subgraphs_rec_check(Graph& g, std::vector<ld>& m, std::bitset<MAXN> cur, std::bitset<MAXN> left, std::bitset<MAXN> neigh, int remDepth, std::string& attr){
        if(remDepth == 0){
            int cnt_cur = cur.count();
            int cnt_left = left.count();
            for(int i=0; i<=cnt_left; i++){
                m[cnt_cur +i] += C[cnt_left][i];
            }
            return;
        }

        std::bitset<MAXN> cand;
        if(cur.count() == 0)
            cand = left;
        else
            cand = left&neigh;

        int cnt = cand.count();
        if(cnt > 0){
            int v = 0;
            while(cand[v] == 0) v++; // v is first form cand
            left.reset(v);
            count_subgraphs_rec_check(g, m, cur, left , neigh, remDepth-1, attr);
            for(auto u : g.adj[v]){
                neigh.set(u.first);
            }
            cur.set(v);
            count_subgraphs_rec_check(g, m, cur, left, neigh, remDepth-1, attr);
        }
        else{
            for(int i=0; i<g.n; i++){
                if(cur[i] && g.node_attr[i]==attr){
                    m[cur.count()]++;
                    return;
                }
            }
        }
    }

    void count_subgraphs_check(Graph& g, std::vector<ld>& m, int remDepth, std::string attr){
        if(g.n >= MAXN){ // don't recur
    	    for(int k=0; k<(int)m.size(); k++){
        		long double lognchoosek = 0;
        		for(int k2=1; k2<=k; k2++) lognchoosek += log(g.n+1-k2)-log(k2);
        		m[k] = lognchoosek;
    	    }
    	    return;
        }
        std::bitset<MAXN> left(0);
        for(int i=0; i<g.n; i++){
            if(g.node_attr[i] <= attr)
                left.set(i);
        }
        count_subgraphs_rec_check(g, m, std::bitset<MAXN>(0), left, std::bitset<MAXN>(0), remDepth, attr);
        for(int k=0; k<(int)m.size(); k++) m[k] = log(m[k]);
    }



    void update_w_tilde(Graph& tau){
        using namespace std;
        unordered_map<string, int> s_set;
        m.resize(max((int)m.size(), tau.n+1));
        t.resize(max((int)m.size(), tau.n+1));
        for(int i=0; i<tau.node_attr.size(); i++){
            s_set[tau.node_attr[i]]++;
        }
        ld cnt_s = 0;
        for(auto pr : s_set){
            string attr = pr.first;
            int cnt_feas = 0;
            for(int i=0; i<tau.n; i++){
                if(tau.node_attr[i] <= attr)
                    cnt_feas++;
            }
            vector<ld> m_tmp(cnt_feas+1, 0.0);
            count_subgraphs_check(tau, m_tmp, maxdepth, attr);
            for(int i=1; i<=cnt_feas; i++){ // subgraph of i nodes
                t[i][attr]++;
                m[i][attr] = logsumexp(m[i][attr], m_tmp[i]);
            }
        }
    }

    double sample_and_update(int to_sample){
        using namespace std;

        FILE* cfile = fopen(infile.c_str(), "r");
        if(VERBOSE) cerr << "Starting sampling " << to_sample << " transactions" << endl;

        vector<int> new_sampled_idxs;
        while(sampled_idxs.size() < to_sample){
            int s = mt()%n;
            if(m.empty()) s = s%(5*to_sample); //speed up first pass
            sampled_idxs.push_back(s);
            new_sampled_idxs.push_back(s);
        }
        sort(new_sampled_idxs.begin(), new_sampled_idxs.end());
        int cnt = 0, id = 0;

        char array[100];
        char* err = fgets ( array, 100, cfile );
        while ( !feof ( cfile ) && id < new_sampled_idxs.size()) {
            Graph tau;
            if(cnt == new_sampled_idxs[id]) readGraph(cfile, tau, true);
            else skipGraph(cfile); // faster
            while(cnt == new_sampled_idxs[id]){
                update_w_tilde(tau);
                id++;
            }
            cnt++;
        }
        fclose(cfile);
        // update rademacher bound
        vector<pair<ld, ld >> data;
        for(int i=1; i<t.size() && i <= maxs; i++){
            if(VERBOSE > 1) cerr << "Size " << i << endl;
            for(auto pr : t[i]){
                string attr = pr.first;
                if(VERBOSE > 1) cerr << " "<< attr << ": " << t[i][attr] << " " << m[i][attr] << endl;
                data.push_back( {min(log(2.0L)*t[i][attr], m[i][attr]), t[i][attr]/(2.0L*to_sample*to_sample)} );
            }
        }


        nlopt_opt opt_prob = nlopt_create(NLOPT_LN_COBYLA, 1);
    	nlopt_set_min_objective(opt_prob, opt_fun, &data);
        double lb[1] = {1.0};
    	nlopt_set_lower_bounds(opt_prob, lb);
    	const double xtol_abs = 1e-4;
    	nlopt_set_xtol_abs(opt_prob, &xtol_abs);
    	nlopt_set_ftol_abs(opt_prob, 1e-6);
    	// Set initialization point
        double x[1]; x[0] = 200.0;
        double min_value;
        nlopt_result result = nlopt_optimize(opt_prob, x, &min_value);
        if(VERBOSE) cerr << "R_S: "<<min_value << " at s = " << x[0]<< endl;
        return min_value; // return rademacher complexity bound
    }

    int next_sample_size(double eta, double eps, int last){
        return std::min(n, (int)ceil((1.1*eta*eta/(eps*eps))*last)); // give 10% more samples
    }



    int n, maxs, maxdepth;
    std::mt19937 mt;
    std::string infile;
    std::vector<std::unordered_map<std::string, ld>> m;
    std::vector<std::unordered_map<std::string, int>> t;
    std::vector<int> sampled_idxs;
    std::vector<std::vector<ld>> C;
};











struct UnlabeledRademacherComputer{
public:
    UnlabeledRademacherComputer(std::string _file, int _n, int _maxs, int seed, int _maxdepth){
        using namespace std;
        n = _n;
        infile = _file;
        maxdepth = _maxdepth;
        mt = std::mt19937( (seed >= 0 ? seed : time(0) ) );
        if(_maxs < 0 || _maxs > MAX_P_SIZE) maxs = MAX_P_SIZE;
        else maxs = _maxs;

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

        // precompute patterns
        vector<vector<pair<int, string>>> v(3);
        v[0] = {{1, "0"}}; v[1] = {{0, "0"},{2, "0"}}; v[2] = {{1, "0"}};
        Graph G_1 = {3, vector<std::string>(3, "0"), v}; vflib::ARGraph<string, string> G_1a = graph2argraph(G_1);
        v[0] = {{1, "0"},{2,"0"}}; v[1] = {{0, "0"},{2, "0"}}; v[2] = {{0,"0"},{1, "0"}};
        Graph G_2 = {3, vector<std::string>(3, "0"), v}; vflib::ARGraph<string, string> G_2a = graph2argraph(G_2);

        v = vector<vector<pair<int, string>>>(4);
        v[0] = {{1, "0"}}; v[1] = {{0, "0"},{2, "0"}}; v[2] = {{1, "0"},{3,"0"}}; v[3] = {{2,"0"}};
        Graph G_3 = {4, vector<std::string>(4, "0"), v}; vflib::ARGraph<string, string> G_3a = graph2argraph(G_3);
        v[0] = {{3, "0"}}; v[1] = {{3, "0"}}; v[2] = {{3, "0"}}; v[3] = {{0,"0"},{1,"0"},{2,"0"}};
        Graph G_4 = {4, vector<std::string>(4, "0"), v}; vflib::ARGraph<string, string> G_4a = graph2argraph(G_4);
        v[0] = {{1, "0"},{3,"0"}}; v[1] = {{0, "0"},{2, "0"}}; v[2] = {{1, "0"},{3,"0"}}; v[3] = {{2,"0"},{0,"0"}};
        Graph G_5 = {4, vector<std::string>(4, "0"), v}; vflib::ARGraph<string, string> G_5a = graph2argraph(G_5);
        v[0] = {{1, "0"},{2,"0"}}; v[1] = {{0, "0"},{2, "0"}}; v[2] = {{0,"0"},{1, "0"},{3,"0"}}; v[3] = {{2,"0"}};
        Graph G_6 = {4, vector<std::string>(4, "0"), v}; vflib::ARGraph<string, string> G_6a = graph2argraph(G_6);
        v[0] = {{1, "0"},{3,"0"},{2,"0"}}; v[1] = {{0, "0"},{2, "0"}}; v[2] = {{1, "0"},{3,"0"},{0,"0"}}; v[3] = {{2,"0"},{0,"0"}};
        Graph G_7 = {4, vector<std::string>(4, "0"), v}; vflib::ARGraph<string, string> G_7a = graph2argraph(G_7);
        v[0] = {{1, "0"},{3,"0"},{2,"0"}}; v[1] = {{0, "0"},{2, "0"},{3,"0"}}; v[2] = {{1, "0"},{3,"0"},{0,"0"}}; v[3] = {{2,"0"},{0,"0"},{1,"0"}};
        Graph G_8 = {4, vector<std::string>(4, "0"), v}; vflib::ARGraph<string, string> G_8a = graph2argraph(G_8);


        patterns = {vflib::ARGraph<string, string>(), G_1a, G_2a, G_3a, G_4a, G_5a, G_6a, G_7a, G_8a};

    }

    double get_epsilon(double delta){
        double r_tilde = sample_and_update(n);
        double eta = 2.0*r_tilde + sqrt(2.0*log(2.0/delta)/n);
        return eta;
    }



    void progressive_sampling(double eps, double delta){
        int samp_size = std::min(n, (int)(8.0*log(2.0/delta)/(eps*eps)) );
        while(true){
            double r_tilde = sample_and_update(samp_size);
            double eta = 2.0*r_tilde + sqrt(2.0*log(2.0/delta)/samp_size);
            if(VERBOSE) std::cerr << "Eta " << eta << std::endl;
            if(eta < eps || samp_size == n){
                break;
            }
            samp_size = next_sample_size(eta, eps, samp_size);
        }
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


private:

    void count_subgraphs_rec_check(Graph& g, std::vector<ld>& m, std::bitset<MAXN> cur, std::bitset<MAXN> left, std::bitset<MAXN> neigh, int remDepth){
        if(remDepth == 0){
            int cnt_cur = cur.count();
            int cnt_left = left.count();
            for(int i=0; i<=cnt_left; i++){
                m[cnt_cur + i] += C[cnt_left][i];
            }
            return;
        }

        std::bitset<MAXN> cand;
        if(cur.count() == 0)
            cand = left;
        else
            cand = left&neigh;

        int cnt = cand.count();
        if(cnt > 0){
            int v = 0;
            while(cand[v] == 0) v++; // v is first form cand
            left.reset(v);
            count_subgraphs_rec_check(g, m, cur, left , neigh, remDepth-1);
            for(auto u : g.adj[v]){
                neigh.set(u.first);
            }
            cur.set(v);
            count_subgraphs_rec_check(g, m, cur, left, neigh, remDepth-1);
        }
        else{
            m[cur.count()]++;
            return;
        }
    }

    void count_subgraphs_check(Graph& g, std::vector<ld>& m, int remDepth){
        if(g.n >= MAXN){
            for(int k=0; k<(int)m.size(); k++){
                long double lognchoosek = 0;
                for(int k2=1; k2<=k; k2++) lognchoosek += log(g.n+1-k2)-log(k2);
                m[k] = lognchoosek;
            }
            return;
        }
        std::bitset<MAXN> left(0);
        for(int i=0; i<g.n; i++){
            left.set(i);
        }
        count_subgraphs_rec_check(g, m, std::bitset<MAXN>(0), left, std::bitset<MAXN>(0), remDepth);
        for(int k=0; k<(int)m.size(); k++) m[k] = log(m[k]);
    }




    void update_w_tilde(Graph& tau){
        using namespace std;

        unordered_map<string, int> s_set;
        m.resize(max((int)m.size(), tau.n+1), vector<ld>(30));
        t.resize(max((int)m.size(), tau.n+1), vector<int>(30));

        for(int j=1; j<9; j++){
            Matcher ma = Matcher(tau, patterns[j]);
            if(ma.isSubgraph()){
                vector<ld> m_tmp(tau.n+1);
                count_subgraphs_check(tau, m_tmp, maxdepth);
                for(int sz=1; sz<=tau.n; sz++){ // subgraph of i nodes
                    t[sz][j]++;
                    m[sz][j] = logsumexp(m[sz][j], m_tmp[sz]); //C[tau.n][sz];
                    if(sz < 16) m[sz][j] = min(m[sz][j], log(n_connected_graphs[sz]));
                }
            }
        }
    }

    double sample_and_update(int to_sample){
        using namespace std;
        FILE* cfile = fopen(infile.c_str(), "r");

        if(VERBOSE) cerr << "Starting sampling " << to_sample << " transactions" << endl;
        vector<int> new_sampled_idxs;
        while(sampled_idxs.size() < to_sample){
            int s = mt()%n;
            if(m.empty()) s = s%(5*to_sample); //speed up first pass
            sampled_idxs.push_back(s);
            new_sampled_idxs.push_back(s);
        }
        sort(new_sampled_idxs.begin(), new_sampled_idxs.end());

        int cnt = 0, id = 0;

        char array[100];
        char* err = fgets ( array, 100, cfile );
        while ( !feof ( cfile ) && id < new_sampled_idxs.size()) {
            Graph tau;
            if(cnt == new_sampled_idxs[id]) readGraph(cfile, tau, true);
            else skipGraph(cfile); // faster

            while(cnt == new_sampled_idxs[id]){
                update_w_tilde(tau);
                id++;
            }
            cnt++;
        }
        fclose(cfile);

        // update rademacher bound
        vector<pair<ld, ld >> data;
        for(int sz=0; sz<t.size() && sz <= maxs; sz++){
            if(VERBOSE > 1) cerr << "Size " << sz << endl;
            for(int i=1; i<3; i++){
                data.push_back( {1.0L, t[sz][i]/(2.0L*to_sample*to_sample)} );
            }
            for(int i=3; i<9; i++){
                if(VERBOSE > 1) cerr << " "<< i << ": " << t[sz][i] << " " << m[sz][i] << endl;
                data.push_back( {min(log(2.0L)*t[sz][i], m[sz][i]), t[sz][i]/(2.0L*to_sample*to_sample)} );
            }
        }


        nlopt_opt opt_prob = nlopt_create(NLOPT_LN_COBYLA, 1);
    	nlopt_set_min_objective(opt_prob, opt_fun, &data);
        double lb[1] = {1.0};
    	nlopt_set_lower_bounds(opt_prob, lb);
    	const double xtol_abs = 1e-4;
    	nlopt_set_xtol_abs(opt_prob, &xtol_abs);
    	nlopt_set_ftol_abs(opt_prob, 1e-6);
    	// Set initialization point
        double x[1]; x[0] = 200.0;
        double min_value;
        nlopt_result result = nlopt_optimize(opt_prob, x, &min_value);
        if(VERBOSE) cerr << "R_S: "<<min_value << " at s = " << x[0]<< endl;
        return min_value; // return rademacher complexity bound
    }

    int next_sample_size(double eta, double eps, int last){
        return std::min(n, (int)ceil((1.1*eta*eta/(eps*eps))*last));
    }


    int n, maxs, maxdepth;
    std::mt19937 mt;
    std::string infile;
    std::vector<vflib::ARGraph<std::string, std::string>> patterns;
    std::vector<std::vector<ld>> m;
    std::vector<std::vector<int>> t;
    std::vector<int> sampled_idxs;
    std::vector<std::vector<ld>> C;
    const ld n_connected_graphs[16] = {1.0, 1.0, 1.0, 2.0, 6.0, 21.0, 112.0, 853.0, 11117.0, 261080.0, 11716571.0, 1006700565.0, 164059830476.0, 50335907869219.0, 29003487462848061.0, 31397381142761241960.0};
};



#endif
