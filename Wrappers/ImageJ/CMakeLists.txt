cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
PROJECT (CCPiQuantificationImageJPlugin)

FIND_PACKAGE(JNI REQUIRED)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

INCLUDE_DIRECTORIES(../Core)
INCLUDE_DIRECTORIES(${JNI_INCLUDE_DIRS})

ADD_LIBRARY(CCPiQuantificationImageJPlugin_Cpp SHARED CCPiAccessibleVolume_wrap.cxx)

FIND_PACKAGE(Java REQUIRED)
INCLUDE(UseJava)


IF(NOT ImageJ)
FIND_JAR(ImageJ ij.jar)
ENDIF(NOT ImageJ)

set(CMAKE_JAVA_INCLUDE_PATH ${ImageJ})

ADD_JAR(CCPiQuantificationImageJPlugin_Java 	AccessibleVolumePluginFilter_.java
							CCPiAccessibleVolume.java
							CCPiAccessibleVolumeInputImages.java
							CCPiAccessibleVolumeITKImpl.java
							CCPiAccessibleVolumeJNI.java
							CCPiConsoleUserInterface.java
							CCPiImageJUserInterface.java
							CCPiUserApplicationInterface.java
							MapType.java
							SimpleHistogramThresholding_.java
							FloatVector.java
							CCPiSimpleHistogramThresholdingITKImplUnsignedChar.java
							LabelQuantification_.java
							CCPiLabelQuantificationITKImplUnsignedChar.java
							CCPiLabelQuantificationResult.java
							StringVector.java
							IntList.java							
		OUTPUT_NAME		CCPiQuantificationImageJ_Plugin) 						 
INSTALL_JAR(CCPiQuantificationImageJPlugin_Java plugins/CCPi)
include_directories(${CCPiCore_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(CCPiQuantificationImageJPlugin_Cpp CCPiCore ${CCPI_ITK_LIBS} ${CCPI_VTK_LIBS})							
INSTALL(TARGETS CCPiQuantificationImageJPlugin_Cpp
		DESTINATION .)
