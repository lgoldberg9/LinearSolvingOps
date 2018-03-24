all: extrema

extrema: extrema.c
	gcc -g -Wall -std=c11 -o extrema extrema.c golden_section_search.c -lm

clean:
	rm ./*~ extrema
