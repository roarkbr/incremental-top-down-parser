include ../Makedefs

LIBINSTALL=$(INSTDIR)/lib
SRC= $(SRCDIR)/lib 
SRCS=parser.c putil.c io-util.c
OBJS=parser.o putil.o io-util.o

all: $(PARLIB)

$(PARLIB): $(OBJS)
	ar rc $@ $(OBJS)

$(OBJS): $(INCS) Makefile ../Makedefs
