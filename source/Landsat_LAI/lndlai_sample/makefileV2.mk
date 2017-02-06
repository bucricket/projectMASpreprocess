#############################################################################
# !make
#
# makefile name: makefile (for ledaps)
#
##!END
#############################################################################
CC = gcc

TARGET = lndlai_sample
BASE = /Users/mschull/umdGD/pyDisALEXI
INCDIR=$(BASE)/processData/source/include/
LIBDIR=$(BASE)/processData/source/lib/

INC =I. -I$(INCDIR)
#LIB =-L$(LIBDIR)

INCS = -I$(INC)

LIBS = -L$(LIBDIR) -lhdfeos -lGctp -lmfhdf -ldf -lz -lsz -ljpeg -lm \

OBJ = lndlai_sample.o lndlai_util.o

all: $(TARGET)

$(OBJ) : $(INC_FILES)

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIBS)  -o $(TARGET)

#
# Rules
#
.c.o: $(INC_FILES)
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(INCS) -c $< -o $@

#******************* End of make file *******************************


