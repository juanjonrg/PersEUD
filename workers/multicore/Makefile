CC=gcc
FLAGS=-O3 -fopenmp
LIBS=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lzmq -lczmq -lcjson

gradient_multicore: gradient_multicore_v1.cpp
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f gradient_multicore
