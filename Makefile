CFLAGS+=-g -ggdb3
CFLAGS+= -O5
LDFLAGS=-lm
PYTHON=python			#Name of Python executable

all: gauss_solve libgauss.so

OBJS = gauss_solve.o main.o helpers.o
gauss_solve.o : gauss_solve.h
helpers.o: helpers.h

gauss_solve : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS)

libgauss.so: $(OBJS)
	$(CC) -shared $(OBJS) -I/usr/include/python3.12 -o $@ -fPIC

check: check_gauss_solve check_ctype_wrapper

check_gauss_solve: gauss_solve
	./$<

check_ctype_wrapper: gauss_solve.py libgauss.so
	$(PYTHON) ./$<

clean: FORCE
	@-rm gauss_solve *.o
	@-rm *.so

FORCE:
