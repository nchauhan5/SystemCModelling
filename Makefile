# Makefile for ECPS 203, Assignment 5 (by Zhongqi Cheng)
# rename the file to Makefile if you want to use it

SYSTEMC = /opt/pkg/systemc-2.3.1

INCLUDE = -I. -I$(SYSTEMC)/include
LIBRARY = $(SYSTEMC)/lib-linux64
CFLAG = $(INCLUDE) -L$(LIBRARY) -Xlinker -R -Xlinker $(LIBRARY) -lsystemc -O2

CC = g++
RM = rm -f

TARGET	= Canny

all: $(TARGET)

test: Canny
	./Canny

clean:
	$(RM) *.o $(TARGET) 

Canny: Canny.cpp
	$(CC) $^ $(CFLAG) -o $@
