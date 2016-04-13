LIBS=-L$(CURDIR) -L/usr/local/lib -L/usr/local/lib -L../flint2/ -L../m4ri/.libs/ -lm4ri -lflint -lmpfr -lgmp -lm -lpthread -Wl,-rpath,/home/best/m4ri/.libs/ -Wl,-rpath,/home/best/flint2/
INCS=-I$(CURDIR) -I/usr/local/include -I/usr/local/include -I../flint2/ -I../m4ri/m4ri/ -I../m4ri/
CFLAGS=-pedantic -Wall -O2 -funroll-loops --std=c99

all:
	$(CC) $(INCS) $(CFLAGS) bw.c $(LIBS) 
