
gccopt= gcc -O2 -funroll-loops
gccdeb= gcc -ggdb -Wall -pedantic
gppdeb = g++ $(CCDEB) 
#CCC = $(gccopt)
CCC = $(gccdeb)

CFLAGS=
LINUXLIBS= -lm -lpthread

LIBS= $(LINUXLIBS) 

CCCLNFLAGS = $(LIBS)

CCCFLAGS = 

PROG = rpower
PROG_OBJ = bin/mainrpower.o bin/power.o


all: bin/$(PROG)

bin/%.o: src/%.c
	@echo compiling $*.c with $(CCC) $(CCCFLAGS)
	@$(CCC) $(CCCFLAGS) -c $< -o $@

bin/$(PROG): $(PROG_OBJ)
	$(CCC) $(CCCFLAGS) -o bin/$(PROG) $(PROG_OBJ) $(CCCLNFLAGS)

clean:
	rm bin/*
