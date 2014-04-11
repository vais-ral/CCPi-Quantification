# Makefile for VolView plugins
#
# David Worth STFC July 2012

CC=g++

VOLVIEW=/home/djw/CCPi/VolViewSrc

VV_INCLUDE=-I$(VOLVIEW)/VolViewLibraries/KWVolView/Plugins

ITK_INCLUDE=-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Code/Common \
-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Code/BasicFilters/ \
-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Code/IO \
-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Code/Review \
-I$(VOLVIEW)/VolViewPlatformComplete/build/ITK

VXL_INCLUDE=-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Utilities/vxl/core \
-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Utilities/vxl/vcl \
-I$(VOLVIEW)/VolViewPlatformComplete/build/ITK/Utilities/vxl/vcl \
-I$(VOLVIEW)/VolViewPlatformComplete/build/ITK/Utilities/vxl/core

VTK_INCLUDE=-I$(VOLVIEW)/VolViewExternalLibraries/VTK/Graphics \
-I$(VOLVIEW)/VolViewExternalLibraries/VTK/Common \
-I$(VOLVIEW)/VolViewExternalLibraries/VTK/Filtering \
-I$(VOLVIEW)/VolViewExternalLibraries/VTK/Imaging \
-I$(VOLVIEW)/VolViewExternalLibraries/VTK/IO \
-I$(VOLVIEW)/VolViewPlatformComplete/build/VTK

APP_INCLUDE=-I$(VOLVIEW)/VolViewExternalLibraries/VTK/IO/ \
$(VTK_INCLUDE) $(ITK_INCLUDE) $(VXL_INCLUDE) $(VV_INCLUDE)


CFLAGS=$(VV_INCLUDE) $(ITK_INCLUDE) $(VXL_INCLUDE) $(VTK_INCLUDE) -fPIC -Wall -std=c++98 -g
# -fopenmp
# Flags for compiler optimisation
#CFLAGS=$(VV_INCLUDE) $(ITK_INCLUDE) $(VXL_INCLUDE) $(VTK_INCLUDE) -fPIC -Wall -std=c++98 -O3 -#D_FORTIFY_SOURCE=0

LIB_DIR=$(VOLVIEW)/VolViewPlatformComplete/build/bin

LDFLAGS=-L$(LIB_DIR) $(LIB_DIR)/libITKIO.so.3.2.0 \
$(LIB_DIR)/libITKCommon.so.3.2.0 $(LIB_DIR)/libITKBasicFilters.so.3.2.0 \
$(LIB_DIR)/libITKAlgorithms.so.3.2.0 $(LIB_DIR)/libITKNrrdIO.so.3.2.0 \
$(LIB_DIR)/libitkjpeg12.so.3.2.0 $(LIB_DIR)/libitkjpeg16.so.3.2.0 \
$(LIB_DIR)/libgdcm.so $(LIB_DIR)/libgdcmjpeg8.so $(LIB_DIR)/libgdcmjpeg12.so \
$(LIB_DIR)/libgdcmjpeg16.so $(LIB_DIR)/libvtkopenjpeg.so.1.3.0 \
$(LIB_DIR)/libvtkopenjp3dvm.so $(LIB_DIR)/libitkpng.so.3.2.0 \
$(LIB_DIR)/libitktiff.so.3.2.0 $(LIB_DIR)/libitkjpeg8.so.3.2.0 \
$(LIB_DIR)/libITKSpatialObject.so.3.2.0 $(LIB_DIR)/libITKMetaIO.so.3.2.0 \
$(LIB_DIR)/libitkopenjpeg.so.1.3.0 $(LIB_DIR)/libitkopenjp3dvm.so \
$(LIB_DIR)/libITKDICOMParser.so.3.2.0 $(LIB_DIR)/libITKEXPAT.so.3.2.0 \
$(LIB_DIR)/libITKniftiio.so.3.2.0 $(LIB_DIR)/libITKznz.so.3.2.0 \
$(LIB_DIR)/libitkzlib.so.3.2.0 $(LIB_DIR)/libITKNumerics.so.3.2.0 \
$(LIB_DIR)/libITKStatistics.so.3.2.0 $(LIB_DIR)/libitkvnl_inst.so.3.2.0 \
$(LIB_DIR)/libitkvnl_algo.so.3.2.0 $(LIB_DIR)/libitkvnl.so.3.2.0 \
$(LIB_DIR)/libitkv3p_netlib.so.3.2.0 $(LIB_DIR)/libitkvcl.so.3.2.0 \
$(LIB_DIR)/libvtkImaging.so $(LIB_DIR)/libvtkFiltering.so \
$(LIB_DIR)/libvtkIO.so  $(LIB_DIR)/libitksys.so.3.2.0 -lpthread -ldl -lm 
# -fopenmp

APP_LDFLAGS=-L$(LIB_DIR)  $(LIB_DIR)/libvtkCommon.so \
$(LIB_DIR)/libvtkGraphics.so $(LIB_DIR)/libvtkImaging.so $(LIB_DIR)/libitkvnl.so \
$(LIB_DIR)/libITKCommon.so 

vvVTKDistanceMap.so : vvVTKDistanceMap.o 
	$(CC) -fPIC -Wall -shared -Wl,--no-undefined -Wl,-soname,$@ -o $@ $^ $(LDFLAGS)
# Use option below in line above to get linker to tell you what symbols are 
# missing rather than find out at run time
# -Wl,--no-undefined

vvVTKDistanceMap.o : vvVTKDistanceMap.cpp
	$(CC) $(CFLAGS) -c $<

clean :
	rm -f *.o *.so*
