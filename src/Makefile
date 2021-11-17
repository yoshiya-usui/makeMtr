CXX           = icpc
CC            = icpc
CXXFLAGS      = -O2 \
                -DNDEBUG 
DEST          = ./
OBJS          = ObservationPoint.o \
                ObservationLine.o \
                main.o
PROGRAM       = makeMtr

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:;		rm -f *.o *~ $(PROGRAM)
