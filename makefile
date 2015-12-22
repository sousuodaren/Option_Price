#CC = /opt/intel/composerxe/bin/icc
CC=gcc -lm
#INC_PATH = -I/opt/intel/composerxe_mic/compiler/include/
INC_PATH = -I/usr/include/


BSDE.out: BSDE.c
	$(CC) BSDE.c -o  BSDE.out  -lpthread
clean:
	rm -f *.o *.out

