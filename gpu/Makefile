OBJS = main.o InitFuncs.o RunSimulation.o OutputFuncs.o

###################################################
#   Compilers    #
###################################################
COMP_GPU=nvcc
COMP=c++
CFLAGS=-O3 -w
FLAGS_GPU=-gencode arch=compute_13,code=sm_13 -gencode arch=compute_20,code=sm_20
NAME=cru3d
XERCES_PATH=../xerces-build
LIBS=xerces-c

##########################################
#               main program             #
##########################################
all:	$(OBJS)
		$(COMP) $(OBJS) $(CFLAGS) -L$(XERCES_PATH)/lib -Wl,-rpath=$(XERCES_PATH)/lib -lxerces-c -L/usr/local/cuda/lib64 -lcudart -o $(NAME)  

clean:
	rm -f *.o $(NAME)

###########################################
#               subroutines               #
###########################################

main.o: main.cpp  parameters.h
	$(COMP_GPU) -c $(CFLAGS) $(INCLUDES) $<

InitFuncs.o: InitFuncs.cu  parameters.h DeviceStructs.h
	$(COMP_GPU) -c $(CFLAGS) -I$(XERCES_PATH)/include $(FLAGS_GPU) $(INCLUDES) $<

RunSimulation.o: RunSimulation.cu  parameters.h DeviceStructs.h
	$(COMP_GPU) -c $(CFLAGS) $(FLAGS_GPU) $(INCLUDES) $< 

OutputFuncs.o: OutputFuncs.cpp  parameters.h
	$(COMP_GPU) -c $(CFLAGS) $(INCLUDES) $<

