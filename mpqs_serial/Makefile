CC = mpicc

# Debugging options
DBG_FLAGS = -Wall -O2 -ggdb

# Performance options
PERF_FLAGS = -Wall -O3 -fomit-frame-pointer

# Profiling options
PROF_FLAGS = -Wall -O2 -ggdb -pg

all:			debug

performance:		
			$(CC) $(PERF_FLAGS) -o qs qs.c -lm -lgmp

debug:			
			$(CC) $(DBG_FLAGS) -o qs qs.c -lm -lgmp

profile:		
			$(CC) $(PROF_FLAGS) -o qs qs.c -lm -lgmp

clean:			
			rm -f core qs *.o *# *~ *.out
