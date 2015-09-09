OBJS = main.o

###################################################
#   Compilers    #
###################################################
COMP=g++
CFLAGS=-O3
NAME=linescan

##########################################
#               main program             #
##########################################
all:	$(OBJS)
	$(COMP) $(OBJS) $(CFLAGS) -L/home/mwalker/models/CRU3D_CPU/xercesc-dist/lib -Wl,-rpath=/home/mwalker/models/CRU3D_CPU/xercesc-dist/lib -lxerces-c -o $(NAME)

clean:
	rm -f *.o $(NAME)

###########################################
#               subroutines               #
###########################################

main.o: main.cpp  parameters.h
	$(COMP) -c $(CFLAGS) -I.. $<
