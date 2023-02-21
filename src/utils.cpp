#include "utils.hpp"

std::string readstring( FILE *input ) {
	char c;
	char word[100];
	int i = 0;
	do{
		c = (char)fgetc( input );
		if(c == ' ' || c == '\n') break;
		else{
			word[i++] = c;
			word[i] = 0;
		}
	} while(c != EOF);
	return (std::string)(word);
}

int readint ( FILE *input ) {
  char car = fgetc ( input );
  while ( car < '0' || car > '9' ) {
    if ( feof ( input ) )
      return -1;
    car = fgetc ( input );
  }
  int n = car - '0';
  car = fgetc ( input );
  while ( car >= '0' && car <= '9' ) {
    n = n * 10 + car - '0';
    car = fgetc ( input );
  }

  return n;
}

void skipToken(FILE* input){
    char c;
    do{
    	c = (char)fgetc( input );
        if(c == ' ' || c == '\n') break;
    } while(c != EOF);
}

void skipLine(FILE* input){
    char c;
    do{
    	c = (char)fgetc( input );
        if(c == '\n') break;
    } while(c != EOF);
}

char readcommand ( FILE *input ) {
  char car = fgetc ( input );
  while ( car < 'a' || car > 'z' ) {
    if ( feof ( input ) )
      return -1;
    car = fgetc ( input );
  }
  return car;
}


void readGraph(FILE* input, Graph& tau, bool labeled){
	char command;
	int nodessize = 0, edgessize = 0;
	command = readcommand ( input );
	while (command == 'v') {
		int id = readint(input);
		std::string lab = readstring(input);
		if(!labeled) lab = "0";
		tau.node_attr.push_back(lab);
		command = readcommand(input);
	}
    tau.n = tau.node_attr.size();
    tau.adj = std::vector<std::vector<std::pair<int, std::string>>>(tau.n, std::vector<std::pair<int, std::string>>(0));
	while (!feof ( input ) && command == 'e'){
		int u = readint(input);
		int v = readint(input);
		std::string lab = readstring(input);
		if(!labeled) lab = "0";
		tau.adj[u].push_back(make_pair(v, lab));
		tau.adj[v].push_back(make_pair(u, lab));
		command = readcommand(input);
	}
}

void skipGraph(FILE* input){
	char command;
    command = readcommand ( input );
    while (command == 'v') {
            skipLine(input);
            command = readcommand(input);
    }
    while (!feof ( input ) && command == 'e'){
		skipLine(input);
		command = readcommand(input);
    }
}

std::vector<Graph> readInput(FILE* input, bool labeled){
	std::vector<Graph> out;
	char array[100];
	char* c;
	c = fgets ( array, 100, input );

	while ( !feof ( input ) ) {
                Graph tau;
		readGraph( input, tau, labeled);
                out.push_back(tau);
		c = fgets ( array, 100, input );
	}
        return out;
}



/*
bool nextTransaction_(std::ifstream& file, Graph& tau, bool labeled){
        using namespace std;
        bool read = true;
        string line, tmp, tmp2, tmp3;
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

        if(!read || line[0]=='\n') return false;
        return true;
    }
*/
