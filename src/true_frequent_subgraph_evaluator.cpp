#include <vector>
#include <string>
#include "subgraph_sampling.hpp"

using namespace std;

int main(int argc, char *argv[]){
    cout << "Frequent Subgraph Sampler Utility" << endl;
    char opt;
    int maxsize = -1, seed = -1, isRademacher = 0, isLabeled = 0;
    while ( ( opt = getopt ( argc, argv, "m:s:rl" ) ) != -1 ) {
      switch ( opt ) {
        case 'm': maxsize = atoi ( optarg ); break;
        case 's': seed = atoi( optarg ); break;
        case 'r': isRademacher = 1; break;
        case 'l': isLabeled = 1; break;
      }
    }

    if ( argc - optind < 3 || argc - optind > 3 ) {
      cout << "Parameters: [-m size] [-s seed] [-r] [-l] delta n input" << endl;
      return 1;
    }

    double delta = atof(argv[optind+0]);
    int n = atoi(argv[optind+1]);
    string in = argv[optind+2];

    double t1, t2;
    t1 = second();


    cout << "Get max deviation of sample from a distribution (with " << n << " samples)" << endl;
    if(isLabeled){
        if(!isRademacher){
            VCSampler vc(in, isLabeled, seed);
            double err = vc.get_epsilon_evc(delta, maxsize);
            cout << "Eps: " << err << endl;
        }
        else{
            LabeledRademacherComputer rc(in, n, maxsize, seed, 0);
            double err = rc.get_epsilon(delta);
            cout << "Eps: " << err << endl;
        }
    }
    else{
        if(!isRademacher){
            VCSampler vc(in, isLabeled, seed);
            double err = vc.get_epsilon_evc(delta, maxsize);
            cout << "Eps: " << err << endl;
        }
        else{
            UnlabeledRademacherComputer rc(in, n, maxsize, seed, 0);
            double err = rc.get_epsilon(delta);
            cout << "Eps: " << err << endl;
        }
    }


    cout << "Total: took " << second()-t1 << " seconds" <<endl;

}
