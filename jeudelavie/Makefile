CC=mpicc
PG=
CFLAGS= -O3 -Wall $(PG) -I/usr/X11/include
LDLIBS= -lX11 -L/usr/X11/lib $(PG)

e= gvie gvie_cycle

all: $(e)

gvie: graph.o

clean:
	/bin/rm -f $(e) *.o gmon.out
