#############################################################################
# !make
#
# makefile name: makefile (for linux)
#
##!END
#############################################################################
CC = gcc

TARGET = gapfillLAI.exe

INC = -I. -I$(HDFINC) -I$(HDFEOS_INC)

LIB =     -L$(HDFEOS_LIB) -lhdfeos -lGctp \
          -L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz -lsz -lm

OBJ = gapfillLAI.o

all: $(TARGET)

$(OBJ) : $(INC_FILES)

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIB)  -o $(TARGET)

#
# Rules
#
.c.o: $(INC_FILES)
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(INC) -c $< -o $@

#******************* End of make file *******************************
