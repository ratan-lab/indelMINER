CC = gcc

OPTIMIZATIONS=
CFLAGS=

INCDIR = 
LNKDIR = 

# the default is the production build without any profiling information
PROD ?= 1
PROF ?= 0

ifeq ($(PROD), 1)
	CFLAGS += -W -Wformat -Wimplicit -Wreturn-type -Wall \
		      -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
			  -Wcast-align -Winline -Wnested-externs -Wextra \
			  -std=c99 -D_USE_KNETFILE -DNDEBUG -O3 -funroll-all-loops \
			  -minline-all-stringops -momit-leaf-frame-pointer 
else
	CFLAGS += -W -Wformat -Wimplicit -Wreturn-type -Wall -Wunused-variable \
			  -Wunused-parameter -Wreturn-type -Wswitch -Wcast-align \
			  -Winline -Wnested-externs -Wextra -std=c99 -D_USE_KNETFILE \
			  -g -ggdb -Werror
endif

ifeq ($(PROF), 1)
	CFLAGS += -pg
endif
	
all: binaries

binaries: indelminer

indelminer:asserts.h \
	       constants.h \
 		   errors.h errors.c \
		   memalloc.h memalloc.c \
		   resources.h resources.c \
		   files.h files.c \
		   sequences.h sequences.c \
		   slinklist.h slinklist.c \
		   superfasthash.h superfasthash.c \
		   hashtable.h hashtable.c \
		   shared.h shared.c \
		   bamoperations.h bamoperations.c \
		   readaln.h readaln.c \
		   evidence.h evidence.c \
		   globalalign.h globalalign.c \
		   localalign.h localalign.c \
	       alignment.h alignment.c \
		   graph.h graph.c \
		   variant.h variant.c \
		   indelminer.c
	$(CC) $(CFLAGS) -c errors.c
	$(CC) $(CFLAGS) -c memalloc.c
	$(CC) $(CFLAGS) -c strings.c
	$(CC) $(CFLAGS) -c resources.c
	$(CC) $(CFLAGS) -c files.c
	$(CC) $(CFLAGS) -c sequences.c
	$(CC) $(CFLAGS) -c slinklist.c
	$(CC) $(CFLAGS) -c superfasthash.c
	$(CC) $(CFLAGS) -c hashtable.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) shared.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) bamoperations.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) evidence.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) graph.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) variant.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) readaln.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) globalalign.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) localalign.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) alignment.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LNKDIR) \
	-o indelminer \
	errors.o memalloc.o strings.o resources.o files.o sequences.o slinklist.o \
	superfasthash.o hashtable.o shared.o bamoperations.o readaln.o evidence.o \
	graph.o variant.o globalalign.o localalign.o alignment.o \
	indelminer.c -lbam -lz -lm

.PHONY: clean archive

clean:
	rm *.o
	rm indelminer
