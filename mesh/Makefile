OBJS = main.o

###################################################
#   Compilers    #
###################################################
COMP=g++
CFLAGS=-O3
NAME=mesh
TET_PATH=../tetgen1.4.3
XERCES_PATH=../xerces-build
##########################################
#               main program             #
##########################################
all:	$(OBJS)
	$(COMP) $(CFLAGS) $(OBJS) -L$(TET_PATH) -ltet -L$(XERCES_PATH)/lib -Wl,-rpath=$(XERCES_PATH)/lib -lxerces-c -o $(NAME)

clean:
	rm -f *.o $(NAME)

###########################################
#               subroutines               #
###########################################

main.o: main.cpp  parameters.h
	$(COMP) -c $(CFLAGS) -I$(TET_PATH) -I$(XERCES_PATH)/include $<
