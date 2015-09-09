OBJS = main.o

###################################################
#   Compilers    #
###################################################
COMP=g++
CFLAGS=-O3
NAME=linescan
XERCES_DIR=../xerces-build

##########################################
#               main program             #
##########################################
all:	$(OBJS)
	$(COMP) $(OBJS) $(CFLAGS) -L$(XERCES_DIR)/lib -Wl,-rpath=$(XERCES_DIR)/lib -lxerces-c -o $(NAME)

clean:
	rm -f *.o $(NAME)

###########################################
#               subroutines               #
###########################################

main.o: main.cpp  parameters.h
	$(COMP) -c $(CFLAGS) -I$(XERCES_DIR)/include $<
