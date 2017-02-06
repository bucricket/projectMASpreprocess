#############################################################################
# !make
#
# makefile name: makefile 
#
##!END
#############################################################################
CC = gcc 

TARGET = th_intC2floatK.exe

OBJ = th_intC2floatK.o

all: $(TARGET)

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(OBJ) -lm  -o $(TARGET)

#******************* End of make file *******************************


