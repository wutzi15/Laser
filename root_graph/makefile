CC	= icc
CFLAGS	=    -fast #-DNDEBUG  
LDFLAGS	= -lSpectrum   -L/Users/wutzis_spam/Downloads/boost_1_46_1/stage/lib/  -lboost_filesystem -lboost_system # -lboost_program_options -lboost_thread 
OBJ	= graph.C graph_main.cpp
#OBJ	= parser.cpp /Users/wutzis_spam/Downloads/boost_1_46_0/stage/lib/libboost_program_options.a /Users/wutzis_spam/Downloads/boost_1_46_0/stage/lib/libboost_thread.a
OPTFL	=  #-fast  -parallel -par-report1 -restrict -override_limits
ROOTCFLAGS = `root-config --cflags` 
ROOTLFLAGS =  `root-config --libs`
analyse:$(OBJ)
	$(CC) $(CFLAGS)$(ROOTCFLAGS) $(OPTFL) -c graph.C 
	$(CC) $(CFLAGS)$(ROOTCFLAGS) $(OPTFL) -c graph_main.cpp
	$(CC) $(LDFLAGS) $(ROOTLFLAGS) -o graph graph.o graph_main.o
root:
	$(CC) $(CFLAGS)$(ROOTCFLAGS) $(OPTFL) -c root2.C 
	$(CC) $(LDFLAGS) $(ROOTLFLAGS) -o graph graph.o

debug:$(OBJ)
	$(CC) $(CFLAGS) $(OPTFL) -o analyse -DDEBUG $(OBJ) $(LDFLAGS)


.PHONY: all
all:
	$(CC) $(CFLAGS) $(OPTFL) -o analyse $(OBJ) $(LDFLAGS)

.PHONY: install
install:
	cp analyse /usr/local/bin/analysewavelength