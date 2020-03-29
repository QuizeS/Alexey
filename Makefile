CC=gcc
CFLAGS?=-Wall -Wextra -pedantic
GDB_CFLAGS?=-Wall -Wextra -pedantic -O0 -ggdb
LDFLAGS?=-lm
RM?=rm -v

ME4: ME4.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

ME4-gdb: ME4.c
	$(CC) $(GDB_CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	$(RM) *.o ME4
