# shell scripts to install all programs for disALEXI data fusion processing #

export CTOP=`pwd`
export CBIN=`echo $CTOP| sed 's/source/bin/'`

export GEOTIFF_INC=$CTOP/include
export GEOTIFF_LIB=$CTOP/lib
export HDFINC=$CTOP/include
export HDFEOS_INC=$CTOP/include
export TIFFINC=$CTOP/include
export TIFFLIB=$CTOP/lib
export HDFEOS_LIB=$CTOP/lib
export HDFLIB=$CTOP/lib

# install Landsat LAI programs
echo "=================="
echo "Installing Landsat LAI ..."
cd Landsat_LAI
for program in "GeoTiff2ENVI" "lndlai_compute" "lndlai_sample"; do
    cd $program
    rm $program.exe
    make
    cp $program.exe $CBIN
    cd ..
done

# install Landsat DMS program
echo "=================="
echo "Installing Landsat DMS ..."
cd $CTOP/Landsat_DMS
for program in "combine_models"  "get_samples"  "predict_fineT"  "th_intC2floatK"; do
    rm $program.exe
    make -f $program.mk
    cp $program.exe $CBIN
done 

#install MODIS DMS program
echo "=================="
echo "Installing MODIS DMS ..."
cd $CTOP/MODIS_DMS
for program in "mod_combine_models" "mod_get_samples" "mod_predict_fineT"; do
    rm $program.exe 
    make -f $program.mk
    cp $program.exe $CBIN
done 
for program in "modlst_C2K" "modndvi_F2I"; do
    cc $program.c -o $program.exe -lm
    cp $program.exe $CBIN
done

# install MODIS LAI smoothing and gap-filling program
echo "=================="
echo "Installing MODIS LAI smoothing and gap-filling ..."
cd $CTOP/Smooth_Gapfill_LAI
for program in "gapfillLAI" "timesatLAI"; do
    rm $program.exe
    make -f $program.mk
    cp $program.exe $CBIN
done 

# install STARFM program
echo "=================="
echo "Installing STARFM ..."
cd $CTOP/StarFM
rm StarFM.exe
make
cp StarFM.exe $CBIN

# install other utilities
echo "==================" 
echo "Installing cubist ..."
cd $CTOP/Cubist
rm cubist
make
cp cubist $CBIN/

echo "==================" 
echo "Installing other utilities ..."
cd $CTOP/LDOPE_QA
rm jdoy
cc jdoy.c -o jdoy -lm
cp jdoy $CBIN/

echo "Intsalled to $CBIN ..."
echo "==================" 