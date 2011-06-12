CC	=  /usr/local/bin/g++
CFLAGS	=   -fno-builtin  -O3 -D_REENTRANT -pthread -m64 -I/Users/wutzis_spam/root/root/include -std=C++0x   #-shared-intel -arch x86_64 -pipe 
LDFLAGS	= -lSpectrum   -L/Users/wutzis_spam/Downloads/boost_1_46_1/stage/lib/  -lboost_filesystem -lboost_system -lboost_program_options 
OPTFL	=  #-fast  -parallel -par-report1 -restrict -override_limits
ROOTCFLAGS = `root-config --cflags` 
ROOTLFLAGS =  `root-config --libs`
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