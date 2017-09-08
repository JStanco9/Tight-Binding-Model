IDIR = ../include
ODIR = ./obj
LDIR = ../lib
BINDIR = ../bin
CC = g++ -std=c++11
DEBUG = -g
CFLAGS = -c -I$(IDIR) $(DEBUG)

_OBJS = tightBinding.o dataInput.o lattice.o tb_help.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))
_DEPS = tightBinding.hpp dataInput.hpp lattice.hpp tb_help.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LFLAGS = $(DEBUG)
LIBS = -I /home/downloads/armadillo-7.960.1/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

all: tbBands findWeyl plotGap photoCurrent test 

tbBands : $(ODIR)/tbBands.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)
	
findWeyl : $(ODIR)/findWeyl.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

plotGap : $(ODIR)/plotGap.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)
	
photoCurrent : $(ODIR)/photoCurrent.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

test : $(ODIR)/test_tb.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

coords : $(ODIR)/ucCoords.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

$(ODIR)/plotGap.o : plotGap.cpp $(DEPS)
	@$(CC) -o $@ $< $(CFLAGS)

$(ODIR)/tightBinding.o : tightBinding.cpp $(DEPS)
	@$(CC) -o $@ $< $(CFLAGS)	
	
$(ODIR)/TB_help.o : tb_help.cpp $(IDIR)/tb_help.hpp
	@$(CC) -o $@ $< $(CFLAGS)
	
$(ODIR)/dataInput.o : dataInput.cpp $(IDIR)/dataInput.hpp $(IDIR)/lattice.hpp	
	@$(CC) -o $@ $< $(CFLAGS)
	
$(ODIR)/lattice.o : lattice.cpp $(IDIR)/lattice.hpp	
	@$(CC) -o $@ $< $(CFLAGS)

$(ODIR)/%.o : %.cpp $(DEPS)
	@$(CC) -o $@ $< $(CFLAGS) 
	
.PHONY: clean

clean : 
	@rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
	
