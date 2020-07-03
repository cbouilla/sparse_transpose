CFLAGS = -O3 -g -mavx2 -Wall -Wextra -std=c11 -fopenmp
CPPFLAGS = -O3 -g -mavx2 -Wall -Wextra -fopenmp

LDFLAGS = -fopenmp
LDLIBS = -lm

#USE_MKL = yes  # comment this line otherwise
#USE_TBB = yes  # comment this line otherwise

ifdef USE_MKL
# for intel MKL
CFLAGS += -DHAVE_MKL -m64 -I${MKLROOT}/include
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed

# parallel MKL using TBB
#LDLIBS += -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -ltbb -lstdc++ -lpthread -lm -ldl

#parallel MKL using Gomp
#LDLIBS += -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

#parallel MKL using iomp5
LDLIBS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# sequential MKL
#LDLIBS += -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
endif


# TBB
ifdef USE_TBB
CFLAGS += -DHAVE_TBB
CPPFLAGS += -DHAVE_TBB
LDLIBS += -ltbb -lstdc++ -lm -ldl
endif


all: driver #parallel_scalability

driver: driver.o mmio.o mini_spasm.o classical.o simple_sort.o

#parallel_scalability: parallel_scalability.o mmio.o mini_spasm.o simple_sort.o


.PHONY: clean

clean:
	rm -rf driver *.o
	#rm -rf parallel_scalability
