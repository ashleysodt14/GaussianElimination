CFLAGS = -std=c99 -fsignaling-nans -g -ggdb3 -O5 -fPIC
LDFLAGS = -lm
PYTHON = python  # Name of Python executable

# Target for the full program
all: gauss_solve libgauss.so

# Object files
OBJS = gauss_solve.o main.o helpers.o

# gauss_solve depends on object files
gauss_solve: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS)

# Rule for building libgauss.so
libgauss.so: gauss_solve.o
	$(CC) -shared -I/usr/include/python3.12 -o $@ -fPIC gauss_solve.o

# Compile gauss_solve.o
gauss_solve.o: gauss_solve.c gauss_solve.h
	$(CC) $(CFLAGS) -c gauss_solve.c

# Compile main.o (if needed)
main.o: main.c helpers.h
	$(CC) $(CFLAGS) -c main.c

# Compile helpers.o (if needed)
helpers.o: helpers.c helpers.h
	$(CC) $(CFLAGS) -c helpers.c

# Test rules
check: check_gauss_solve check_ctype_wrapper

check_gauss_solve: gauss_solve
	./$<

check_ctype_wrapper: gauss_solve.py libgauss.so
	$(PYTHON) ./$<

# Clean rule
clean:
	rm -f gauss_solve *.o libgauss.so

# Force target to always execute clean
FORCE:
