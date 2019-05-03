
# Assignment 2 pi

CC = mpicc

ppi: ppi.o 
	$(CC) ppi.o -o ppi 

ppi.o: ppi.c
	$(CC) -O3 -xHOST -c ppi.c  -o ppi.o

pi: pi.o 
	$(CC) pi.o -o pi 

pi.o: pi.c
	$(CC) -O3 -xHOST -c pi.c  -o pi.o

ppi_weak: ppi_weak.o
	$(CC) ppi_weak.o -o ppi_weak 

ppi_weak.o: ppi_weak.c
	$(CC) -O3 -xHOST -c ppi_weak.c  -o ppi_weak.o

ppi_mr: ppi_mr.o
	$(CC) ppi_mr.o -o ppi_mr

ppi_mr.o: ppi_mr.c
	$(CC) -O3 -xHOST -c ppi_mr.c  -o ppi_mr.o

clean:
	rm -rf *.o core ppi pi ppi_weak ppi_mr

