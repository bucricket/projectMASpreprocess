#############################################################################
# !make
#
# makefile name: timesat.mk (for landsc1.nascom.nasa.gov)
#
##!END
#############################################################################
F90 = gfortran -m32 -static

HDFLIB2 = .

TARGET = timesatLAI.exe

INC = -I. -I$(HDFINC) -I$(HDFEOS_INC)

LIB = -L. -L$(HDFLIB)/m32 -lmfhdf -ldf -lz -lsz -ljpeg -lhdfeos -lGctp -lm

OBJ = basis.f90 fitgauss.f90 fitlogistic.f90 fungauss.f90 funlogistic.f90 \
      gauss.f90 linlsq.f90 marquardt.f90 modweight.f90 median.f90 \
      phenology.f90 savgol.f90 season.f90 spike.f90 \
      timesatLAI.f90

all: $(TARGET)

# Make the process
$(TARGET) : $(OBJ)
	$(F90) $(OBJ) $(LIB) -o $(TARGET)

#******************* End of make file *******************************


