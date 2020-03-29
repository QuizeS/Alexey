CC=gcc
CFLAGS?=-Wall -Wextra -pedantic
GDB_CFLAGS?=-Wall -Wextra -pedantic -O0 -ggdb
LDFLAGS?=-lm
RM?=rm -v

dME4: dME4.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

dME4-gdb: dME4.c
	$(CC) $(GDB_CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	$(RM) *.o dME4
