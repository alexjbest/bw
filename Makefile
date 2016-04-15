LIBS=-L$(CURDIR) -L/usr/local/lib -L/usr/local/lib -L../flint2/ -L../m4ri/.libs/ -lm4ri -lflint -lmpfr -lgmp -lm -lpthread -Wl,-rpath,/home/best/m4ri/.libs/ -Wl,-rpath,/home/best/flint2/
INCS=-I$(CURDIR) -I/usr/local/include -I/usr/local/include -I../flint2/ -I../m4ri/m4ri/ -I../m4ri/
CFLAGS=-pedantic -Wall -O3 -funroll-loops --std=c99 -g -msse2 -ftree-vectorizer-verbose=0
DEPS = bw.h
OBJ = bw.o qseive.o test.o

%.o: %.c $(DEPS)
	$(CC) $(INCS) $(CFLAGS) -c -o $@ $< $(LIBS)

test: bw.o test.o
	$(CC) $(INCS) $(CFLAGS) -o $@ $^ $(LIBS)

qseive: bw.o qseive.o
	$(CC) $(INCS) $(CFLAGS) -o $@ $^ $(LIBS)
