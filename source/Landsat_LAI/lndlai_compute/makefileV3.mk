#############################################################################
# !make
#
# makefile name: makefile (for ledaps)
#
##!END
#############################################################################
CC = gcc

INCDIR=/Applications/anaconda/envs/pyDMS3/include/
LIBDIR=/Applications/anaconda/envs/pyDMS3/lib/

TARGET = lndlai_compute.exe

INC = -I. -I$(INCDIR) -I$(INCDIR)

LIB = -L$(LIBDIR) -lhdfeos -lGctp \
      -L$(LIBDIR) -lmfhdf -ldf -lz -lsz -ljpeg -lm \

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


