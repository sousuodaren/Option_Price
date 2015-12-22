#CC = /opt/intel/composerxe/bin/icc
CC=gcc -lm
#INC_PATH = -I/opt/intel/composerxe_mic/compiler/include/
INC_PATH = -I/usr/include/


 
BSDE.out: BSDE.o ThetaScheme.o 
	$(CC) -O2 -o BSDE.out BSDE.o ThetaScheme.o 

BSDE.o: BSDE.c ThetaScheme.h
	$(CC) -O2 -c -o BSDE.o BSDE.c 
	
ThetaScheme.o:ThetaScheme.c ThetaScheme.h
	$(CC) -O2 -c -o ThetaScheme.o ThetaScheme.c

clean:
	rm -f *.o *.out

