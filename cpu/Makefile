OBJS = main.o InitFuncs.o RunSimulation.o OutputFuncs.o

###################################################
#   Compilers    #
###################################################
COMP=g++
COMP_MPI=mpic++
CFLAGS=-O3
NAME=cru3d
XERCES_PATH=../xerces-build

##########################################
#               main program             #
##########################################
all:	$(OBJS)
		$(COMP_MPI) $(OBJS) $(CFLAGS) -L$(XERCES_PATH)/lib -Wl,-rpath=$(XERCES_PATH)/lib -lxerces-c -o $(NAME)  

clean:
	rm -f *.o $(NAME)

###########################################
#               subroutines               #
###########################################

main.o: main.cpp  parameters.h
	$(COMP_MPI) -c $(CFLAGS) $(INCLUDES) $<

InitFuncs.o: InitFuncs.cpp  parameters.h
	$(COMP_MPI) -c $(CFLAGS) -I$(XERCES_PATH)/include $(INCLUDES) $<

RunSimulation.o: RunSimulation.cpp  parameters.h
	$(COMP) -c $(CFLAGS) $(INCLUDES) $< 

OutputFuncs.o: OutputFuncs.cpp  parameters.h
	$(COMP_MPI) -c $(CFLAGS) $(INCLUDES) $<

