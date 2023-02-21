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

    if ( argc - optind < 5 || argc - optind > 5 ) {
      cout << "Parameters: [-m size] [-s seed] [-r] [-l] eps delta n input output" << endl;
      return 1;
    }

    double eps = atof(argv[optind]);
    double delta = atof(argv[optind+1]);
    int n = atoi(argv[optind+2]);
    string in = argv[optind+3];
    string out = argv[optind+4];
    double t1, t2;
    t1 = second();

	cout << "Sampling dataset ";
    if(!isRademacher){ // use VC
        cout << "using VC" << endl;
        VCSampler vc(in, isLabeled, seed);
        vc.compute_vc_and_sample(eps, delta, maxsize);
        t2 = second();
        cout << "Bound: took " << t2-t1 << " seconds" <<endl;
        int samps = vc.output_sampled_dataset(out);
        cout << "Sampled " << samps << " transactions" << endl;
        cout << "Sample: took " << second()-t2 << " seconds" <<endl;
    }
    else{ // use rademacher progressive sampling
        cout << "using Rademacher" << endl;
        if(isLabeled){
            LabeledRademacherComputer rc(in, n, maxsize, seed, 0);
            rc.progressive_sampling(eps, delta);
            t2 = second();
            cout << "Bound: took " << t2-t1 << " seconds" <<endl;
            int samps = rc.output_sampled_dataset(out);
            cout << "Sampled " << samps << " transactions" << endl;
            cout << "Sample: took " << second()-t2 << " seconds" <<endl;
        }
        else{
            UnlabeledRademacherComputer rc(in, n, maxsize, seed, 0);
            rc.progressive_sampling(eps, delta);
            t2 = second();
            cout << "Bound: took " << t2-t1 << " seconds" <<endl;
            int samps = rc.output_sampled_dataset(out);
            cout << "Sampled " << samps << " transactions" << endl;
            cout << "Sample: took " << second()-t2 << " seconds" <<endl;
        }
    }



    cout << "Total: took " << second()-t1 << " seconds" <<endl;

}
