#ifndef VF3_WRAPPERS_H
#define VF3_WRAPPERS_H

#include "VFLib.h"
#include "graph.hpp"
#include "VF3SubState.hpp"

#define state_t vflib::VF3SubState<std::string, std::string, std::string, std::string>


static inline vflib::ARGraph<std::string, std::string> graph2argraph(Graph& g){
    vflib::ARGLoader<std::string, std::string>* loader = new vflib::FastStreamARGLoader<std::string, std::string>(g);
    vflib::ARGraph<std::string, std::string> arg(loader);
    delete loader;
    return arg;
}

class Matcher {
public:
    Matcher(vflib::ARGraph<std::string, std::string> g, vflib::ARGraph<std::string, std::string> p){
        targ_graph = g;
        patt_graph = p;
    }

    Matcher(Graph& g, Graph& p){
        vflib::ARGLoader<std::string, std::string>* loader = new vflib::FastStreamARGLoader<std::string, std::string>(g);
        targ_graph = vflib::ARGraph<std::string, std::string> (loader);
        delete loader;

        vflib::ARGLoader<std::string, std::string>* loader2 = new vflib::FastStreamARGLoader<std::string, std::string>(p);
        patt_graph = vflib::ARGraph<std::string, std::string>(loader2);
        delete loader2;
    }

    Matcher(Graph& g, vflib::ARGraph<std::string, std::string> p){
        vflib::ARGLoader<std::string, std::string>* loader = new vflib::FastStreamARGLoader<std::string, std::string>(g);
        targ_graph = vflib::ARGraph<std::string, std::string> (loader);
        delete loader;

        patt_graph = p;
    }

    bool isSubgraph(){
        vflib::MatchingEngine<state_t>* me = new vflib::MatchingEngine<state_t>(true);
        vflib::FastCheck<std::string, std::string, std::string, std::string> check(&patt_graph, &targ_graph);

        if(check.CheckSubgraphIsomorphism())
    	{
    		vflib::NodeClassifier<std::string, std::string> classifier(&targ_graph);
    		vflib::NodeClassifier<std::string, std::string> classifier2(&patt_graph, classifier);
            std::vector<uint32_t> class_patt;
	        std::vector<uint32_t> class_targ;
    		class_patt = classifier2.GetClasses();
    		class_targ = classifier.GetClasses();
    		uint32_t classes_count = classifier.CountClasses();

            vflib::VF3NodeSorter<std::string, std::string, vflib::SubIsoNodeProbability<std::string, std::string>> sorter(&targ_graph);
			std::vector<vflib::nodeID_t> sorted = sorter.SortNodes(&patt_graph);

			state_t s0(&patt_graph, &targ_graph, class_patt.data(), class_targ.data(), classes_count, sorted.data());
			//me->FindAllMatchings(s0);
            me->FindFirstMatching(s0);
            bool ans = ( me->GetSolutionsCount() > 0 );
            delete me;
            return ans;
        }
        else{
            delete me;
            return false;
        }
    }

    int countSubgraphs(){
        return -1;
    }


    bool isIsomorphic(){
        vflib::MatchingEngine<state_t>* me = new vflib::MatchingEngine<state_t>(false);
        vflib::FastCheck<std::string, std::string, std::string, std::string> check(&patt_graph, &targ_graph);

        if(check.CheckIsomorphism())
    	{
    		vflib::NodeClassifier<std::string, std::string> classifier(&targ_graph);
    		vflib::NodeClassifier<std::string, std::string> classifier2(&patt_graph, classifier);
            std::vector<uint32_t> class_patt;
	        std::vector<uint32_t> class_targ;
    		class_patt = classifier2.GetClasses();
    		class_targ = classifier.GetClasses();
    		uint32_t classes_count = classifier.CountClasses();

            vflib::VF3NodeSorter<std::string, std::string, vflib::SubIsoNodeProbability<std::string, std::string>> sorter(&targ_graph);
			std::vector<vflib::nodeID_t> sorted = sorter.SortNodes(&patt_graph);

			state_t s0(&patt_graph, &targ_graph, class_patt.data(), class_targ.data(), classes_count, sorted.data());
			bool ans = me->FindFirstMatching(s0);
            delete me;
            return ans;
        }
        else{
            delete me;
            return false;
        }
    }


private:
    vflib::ARGraph<std::string, std::string> targ_graph, patt_graph;
};




#endif
