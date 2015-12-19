TARGET  = Reduce_knot20
SUBTREE =

SRCs= main.c \
      r250.c \
      io.c \
      random_walk.c \
      knot.c \
      determine_knot_size.c \
      taylor.c \
      alexander_polynomial.c \
      density_prof.c \
      analysis.c \
      numerical.c
      

## OBJs:=$(subst .c,.o,$(SRCs))

OBJs= main.o \
      r250.o \
      io.o \
      random_walk.o \
      knot.o \
      determine_knot_size.o \
      taylor.o \
      alexander_polynomial.o \
      density_prof.o \
      analysis.o \
      numerical.o
 
HDRs= element.h \
      r250.h

CC 	   = cc
#CFLAGS     = $(STDCFLAGS) -g
#CC         = icc
#CFLAGS     = $(STDCFLAGS) -O3 -xW -tpp7 -ipo -ipo_obj
CFLAGS     = $(STDCFLAGS) -Ofast 
#CFLAGS     = -g
LDFLAGS    = -lm 

all: python mc 

say	:
	echo $(SRCs)
	echo $(OBJs)

## Here is how to create the objects:

%.o	: %.c $(HDRs)
	$(CC) $(CFLAGS) -c $< `python-config --includes`

mc	:$(OBJs)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJs) det.o $(LDFLAGS)  `python-config --libs`
python 	: 
	cython det.pyx
	gcc -g -O2 -fpic -c det.c -o det.o `python-config --includes`

clean	:
	rm -f $(OBJs) $(TARGET)
