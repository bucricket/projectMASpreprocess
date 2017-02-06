#############################################################################
# !make
#
# makefile name: makefile (for linux)
#
##!END
#############################################################################
CC = gcc

BASE = /Users/mschull/umdGD/pyDisALEXI
INCDIR=$(BASE)/processData/source/include/
LIBDIR=$(BASE)/processData/source/lib/
#INCDIR=/Applications/anaconda/envs/pyDMS3/include/
#LIBDIR=/Applications/anaconda/envs/pyDMS3/lib/

#CFLAGS=-Wl,-Bstatic

TARGET = GeoTiff2ENVI

INC = -I. -I$(INCDIR)

LIB =  -L$(LIBDIR) -lgeotiff  -ltiff -ljpeg -lz -lm \


OBJ = GeoTiff2ENVI.o

all: $(TARGET)

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(OBJ) $(CFLAGS) $(LIB) -o $(TARGET)

#
# Rules
#
.c.o: $(INC_FILES)
	$(CC) $(INC) $(CFLAGS) -c $< -o $@

#******************* End of make file *******************************
