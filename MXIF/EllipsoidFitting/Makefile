# Makefile for VolView plugins
#
# David Worth STFC May 2012

CC=g++

VOLVIEW=/home/djw/CCPi/VolViewSrc
BUILD_DIR=$(VOLVIEW)/build

VV_INCLUDE=-I$(VOLVIEW)/VolViewLibraries/KWVolView/Plugins

ITK_INCLUDE=-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Code/Common \
-I$(BUILD_DIR)/ITK

VXL_INCLUDE=-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Utilities/vxl/core \
-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Utilities/vxl/vcl \
-I$(BUILD_DIR)/ITK/Utilities/vxl/vcl \
-I$(BUILD_DIR)/ITK/Utilities/vxl/core

VTK_INCLUDE=-I$(VOLVIEW)/VolViewExternalLibraries/VTK/Graphics \
-I$(VOLVIEW)/VolViewExternalLibraries/VTK/Common \
-I$(VOLVIEW)/VolViewExternalLibraries/VTK/Filtering \
-I$(VOLVIEW)/VolViewExternalLibraries/VTK/Imaging \
-I$(BUILD_DIR)/VTK

APP_INCLUDE=-I$(VOLVIEW)/VolViewExternalLibraries/VTK/IO/ \
-I$(VOLVIEW)/VolViewExternalLibraries/Insight/Code/IO \
$(VTK_INCLUDE) $(ITK_INCLUDE) $(VXL_INCLUDE) $(VV_INCLUDE)


CFLAGS=$(VV_INCLUDE) $(ITK_INCLUDE) $(VXL_INCLUDE) $(VTK_INCLUDE) -fPIC -Wall -std=c++98 -g -fopenmp
# Flags for compiler optimisation
#CFLAGS=$(VV_INCLUDE) $(ITK_INCLUDE) $(VXL_INCLUDE) $(VTK_INCLUDE) -fPIC -Wall -std=c++98 -O3 -D_FORTIFY_SOURCE=0

LIB_DIR=$(BUILD_DIR)/bin

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
$(LIB_DIR)/libitksys.so.3.2.0 -lpthread -ldl -lm

APP_LDFLAGS=-L$(LIB_DIR) $(LIB_DIR)/libvtkFiltering.so $(LIB_DIR)/libvtkCommon.so \
$(LIB_DIR)/libvtkGraphics.so $(LIB_DIR)/libvtkImaging.so $(LIB_DIR)/libitkvnl.so \
$(LIB_DIR)/libITKCommon.so $(LIB_DIR)/libvtkIO.so

vvEllipsoidFitting.so : EllipsoidFitting.o vvEllipsoidFitting.o EllipsoidWorker.o
	$(CC) -fPIC -fopenmp -Wall -shared -Wl,-soname,$@ -o $@ $^ $(LDFLAGS)
# Use option below in line above to get linker to tell you what symbols are 
# missing rather than find out at run time
# -Wl,--no-undefined

vvEllipsoidFitting.o : vvEllipsoidFitting.cpp EllipsoidFitting.hpp EllipsoidWorker.hpp
	$(CC) $(CFLAGS) -c $<

EllipsoidFitting.o : EllipsoidFitting.cpp
	$(CC) $(CFLAGS) -c $<

EllipsoidFitting.cpp : EllipsoidFitting.hpp EllipsoidWorker.hpp

EllipsoidWorker.o : EllipsoidWorker.cpp
	$(CC) $(CFLAGS) -c $<

EllipsoidWorker.cpp : EllipsoidWorker.hpp

app: AppEllipsoidFitting.o EllipsoidFitting.o EllipsoidWorker.o
	$(CC) -fopenmp -o $@ $^ $(APP_LDFLAGS)

AppEllipsoidFitting.o: AppEllipsoidFitting.cpp EllipsoidFitting.hpp EllipsoidWorker.hpp
	$(CC) -Wall -std=c++98 -g -fopenmp $(APP_INCLUDE) -c $<

clean :
	rm -f *.o *.so*
