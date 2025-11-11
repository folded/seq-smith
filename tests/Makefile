CC = gcc
CFLAGS = -c -fPIC -I align
LDFLAGS = -shared

SRCS = align/align.c align/alignment.c align/global_align.c align/glocal_align.c align/local_align.c
OBJS = $(patsubst %.c,%.o,$(SRCS))
TARGET = libalign.so

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
