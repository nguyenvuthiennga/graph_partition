PHONY: partition all

all: partition


partition: graph.c partition.c utilities.c
	gcc -O0 -g -o partition graph.c partition.c utilities.c `pkg-config --libs --cflags gsl`

clean:
	rm -f partition

