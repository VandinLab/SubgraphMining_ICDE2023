GCC=c++-11 $(RPATH)

INC=-I./vf3lib -I./../nlopt-2.7.1/install/include -L./../nlopt-2.7.1/install/lib
LIBS=utils.cpp -lnlopt -lm
FLAGS=-Wno-deprecated -O3

all: frequent true_frequent match sampler

frequent: *.hpp *.h *.cpp Makefile
	$(GCC) $(INC) $(FLAGS) -o frequent_subgraph_sampler frequent_subgraph_sampler.cpp $(LIBS)

true_frequent: *.hpp *.h *.cpp Makefile
	$(GCC) $(INC) $(FLAGS) -o true_frequent_subgraph_evaluator true_frequent_subgraph_evaluator.cpp $(LIBS)

match: *.hpp *.h *.cpp Makefile
	 $(GCC) $(INC) $(FLAGS) -o match match.cpp $(LIBS)

sampler: *.hpp *.h *.cpp Makefile
	$(GCC) $(INC) $(FLAGS) -o sampler sampler.cpp $(LIBS)
