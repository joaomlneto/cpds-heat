NVCC 	= /opt/cuda/4.1/bin/nvcc
CUDAPATH = /opt/cuda/4.1

NVCCFLAGS = -I$(CUDAPATH)/include #--ptxas-options=-v 
LFLAGS = -L$(CUDAPATH)/lib64 -lcuda -lcudart -lm

ALL	= heatCUDA
all: $(ALL)

kernels.o: kernels.cu 
	$(NVCC) -c -g $(NVCCFLAGS) $+ $(LFLAGS) -o $@

heatCUDA: heatCUDA.cu kernels.o
	$(NVCC) -g -O2 $(NVCCFLAGS) $+ $(LFLAGS) -o $@

clean:
	rm -fr $(ALL) *.o *.prv *.pcf *.row *.sym *.mpits set-0

