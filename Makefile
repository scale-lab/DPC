#
# Makefile that builds btest and other helper programs for the CS:APP data project
# 
CC = gcc
CFLAGS = -pthread -O -Wall -std=gnu99
LIBS = -lm -lrt

.PHONY: all
all: DPC

DPC: DPC.c
	$(CC) DPC.c $(CFLAGS) $(LIBS) -o DPC

.PHONY: clean
clean:
	rm -f *.o DPC *~


