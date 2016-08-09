###################################################################
#
# freeman.justin@gmail.com 
#
##################################################################

CC=	gcc

CSRC=	./src/

CFLAGS=	-O3 -g -Wall 
#CFLAGS=	-O3 -g -fPIC -Wall

INC=	-I./include

#LFLAGS= -lnetcdf -lxml2 
LFLAGS=  \
	-lnetcdf

COBJ=	$(CSRC)main.o

OBJ=	$(COBJ) 

EXEC=	./bin/MetalstorM

$(EXEC):$(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJ) $(LFLAGS)

$(COBJ) : %.o : %.c
	$(CC) $(INC) $(CFLAGS) -c $< -o $@

clean:
	rm $(COBJ)
	rm $(EXEC)
