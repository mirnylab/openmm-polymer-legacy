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
CFLAGS     = $(STDCFLAGS) -Ofast  -m32
#CFLAGS     = -g
LDFLAGS    = -lm 

all: mc

say	:
	echo $(SRCs)
	echo $(OBJs)

## Here is how to create the objects:

%.o	: %.c $(HDRs)
	$(CC) $(CFLAGS) -c $<

mc	: $(OBJs)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJs) $(LDFLAGS) 

clean	:
	rm -f $(OBJs) $(TARGET)
