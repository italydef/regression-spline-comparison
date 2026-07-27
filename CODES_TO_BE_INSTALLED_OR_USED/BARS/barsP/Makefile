CC=gcc
CFLAGS=-std=gnu89 -Wno-implicit-function-declaration -Wno-implicit-int

com.o: com.c
	$(CC) $(CFLAGS) -c com.c

ranlibsub.o: ranlibsub.c ranlib.h
	$(CC) $(CFLAGS) -c ranlibsub.c

logspline.o: logspline.c logspline.h
	$(CC) $(CFLAGS) -c logspline.c

barsP_utils.o: barsP_utils.c barsP_utils.h
	$(CC) $(CFLAGS) -c barsP_utils.c

barsP_funcs.o: barsP_funcs.c barsP_funcs1.h barsP_funcs.h
	$(CC) $(CFLAGS) -c barsP_funcs.c

barsP.out: barsP.c barsP_funcs.h barsP_funcs.o barsP_utils.o logspline.o ranlibsub.o com.o
	$(CC) -o barsP.out $(CFLAGS) barsP.c barsP_funcs.o barsP_utils.o logspline.o ranlibsub.o com.o -framework Accelerate
