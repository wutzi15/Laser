CC	=  /usr/local/bin/g++
CFLAGS	=   -fno-builtin  -Ofast -mtune=core2 -D_REENTRANT -pthread -m64 -I/Users/wutzis_spam/root/root/include -std=c++0x -fno-pic  -fwhole-program -fomit-frame-pointer  -ffast-math -ffunction-sections   -fopenmp -static -S  -msse4.1 -fexcess-precision=fast -floop-parallelize-all -flto 
LDFLAGS	= -lSpectrum   -L/Users/wutzis_spam/Downloads/boost_1_46_1/stage/lib/  -lboost_filesystem -lboost_system -lboost_program_options -flto 
OPTFL	=  #-fast  -parallel -par-report1 -restrict -override_limits
ROOTCFLAGS = `root-config --cflags` 
ROOTLFLAGS =  `root-config --libs`
.PHONY: root_graph
root_graph:$(OBJ)
	$(CC) $(CFLAGS)$(ROOTCFLAGS) $(OPTFL) -c graph.cpp
	$(CC) $(CFLAGS)$(ROOTCFLAGS) $(OPTFL) -c graph_main.cpp
	$(CC) $(CFLAGS)$(ROOTCFLAGS) $(OPTFL) -c helper.cpp
	$(CC) $(LDFLAGS) $(ROOTLFLAGS) -o root_graph graph.o graph_main.o helper.o



.PHONY: install
install:
	cp root_graph /usr/local/bin/root_graph

.PHONY: remake
	$(CC) $(CFLAGS)$(ROOTCFLAGS) $(OPTFL) -c graph.cpp
	$(CC) $(CFLAGS)$(ROOTCFLAGS) $(OPTFL) -c graph_main.cpp
	$(CC) $(CFLAGS)$(ROOTCFLAGS) $(OPTFL) -c helper.cpp
	$(CC) $(LDFLAGS) $(ROOTLFLAGS) -o root_graph graph.o graph_main.o helper.o