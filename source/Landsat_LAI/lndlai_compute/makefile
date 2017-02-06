#############################################################################
# !make
#
# makefile name: makefile (for ledaps)
#
##!END
#############################################################################
CC = gcc

TARGET = lndlai_compute.exe

INC = -I. -I$(HDFINC) -I$(HDFEOS_INC)

LIB = -L$(HDFEOS_LIB) -lhdfeos -lGctp \
      -L$(HDFLIB) -lmfhdf -ldf -lz -lsz -ljpeg -lm \

OBJ = lndlai_compute.o lndlai_util.o lndlai_rt.o

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


