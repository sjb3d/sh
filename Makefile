
CC=clang
CCFLAGS=-O2 -Wall -Wextra

all: sh_gen

clean:
	$(RM) sh_gen

sh_gen: sh_gen.c
	$(CC) $(CCFLAGS) -o$@ $<
