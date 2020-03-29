CC=gcc
CFLAGS?=-Wall -Wextra -pedantic
GDB_CFLAGS?=-Wall -Wextra -pedantic -O0 -ggdb
LDFLAGS?=-lm -llua
RM?=rm -v

ME4: ME4.c
	$(CC) $(CFLAGS) $(_DBG) $(LDFLAGS) -o $@ $<

ME4-gdb: ME4.c
	$(CC) $(GDB_CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	$(RM) *.o ME4
