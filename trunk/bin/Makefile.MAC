include ../Makedefs

SRCS=tdptrain tdparse
OBJS=tdptrain.o tdparse.o
PROGS=tdptrain.c tdparse.c
SRC= $(SRCDIR)/bin

all: $(SRCS)

tdptrain: tdptrain.o $(INCS) $(LIBS)
	$(CC) -o tdptrain tdptrain.o $(LIBS)

tdparse: tdparse.o $(INCS) $(LIBS)
	$(CC) -o tdparse tdparse.o $(LIBS)

$(PARLIB): $(LIBDEP) $(INCS)
	cd $(INSTDIR)/lib; $(MAKE)

$(OBJS): $(INCS) Makefile ../Makedefs


