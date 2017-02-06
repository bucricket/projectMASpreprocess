#############################################################################
# !make
#
# makefile name: makefile 
#
##!END
#############################################################################
CC = gcc 
#CFLAGS = -ansi -Wall -v -da -Q
CFLAGS = -Wall

TARGET = mod_predict_fineT.exe

INC = -I. -I$(TIFFINC) -I$(GEOTIFF_INC)
#INC = -I. -I/usr/include -I/usr/include/libgeotiff

LIB =     -L$(GEOTIFF_LIB) -lgeotiff \
          -L$(TIFFLIB) -ltiff -lm -ljpeg -lz -lsz\

OBJ = predict_fineT.o sensor.o utility.o

all: $(TARGET)

$(OBJ) : modis.h

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIB)  -o $(TARGET)

#
# Rules
#
.c.o: $(INC_FILES)
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(INC) -c $< -o $@

#******************* End of make file *******************************


