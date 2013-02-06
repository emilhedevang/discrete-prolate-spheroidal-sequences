
dpss.o: dpss.c
	gcc -std=c11 -O3 -Wall -o dpss.o -c dpss.c

dpss: dpss.o
	gcc -o dpss -lhdf5 -lhdf5_hl -llapacke dpss.o
