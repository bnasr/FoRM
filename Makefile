CC=g++
DEPS = model.h auxfunc.h svfpl.h sunearth.h prepost.h calculate.h

PROJNAME=form
PROJNAMESERIAL=formserial

LIBS=-lm

CFLAGS=-O3 -Wall
OMPFLAGS=-fopenmp 

ODIR=obj
ODIRSERIAL=objserial

_OBJ= model.o auxfunc.o svfpl.o sunearth.o prepost.o calculate.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
OBJSERIAL = $(patsubst %,$(ODIRSERIAL)/%,$(_OBJ))

all: $(PROJNAME) $(PROJNAMESERIAL)
	
run: $(PROJNAME)
	./$(PROJNAME) test

runserial: $(PROJNAMESERIAL)
	./$(PROJNAMESERIAL) test

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS)

$(ODIRSERIAL)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
	
$(PROJNAME): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(OMPFLAGS) $(LIBS)

$(PROJNAMESERIAL): $(OBJSERIAL)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o ./$(PROJNAME) $(ODIRSERIAL)/*.o ./$(PROJNAMESERIAL)
