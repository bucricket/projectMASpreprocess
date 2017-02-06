#############################################################################
# !make
#
# makefile name: makefile (for linux)
#
##!END
#############################################################################
CC = cc

TARGET = GeoTiff2ENVI.exe

INC = -I. -I$(GEOTIFF_INC)

LIB = -L. -L$(GEOTIFF_LIB) -lgeotiff  -ltiff -ljpeg -lz -lm


OBJ = GeoTiff2ENVI.o

all: $(TARGET)

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(OBJ) $(LIB)  -o $(TARGET)

#
# Rules
#
.c.o: $(INC_FILES)
	$(CC) $(INC) -c $< -o $@

#******************* End of make file *******************************


