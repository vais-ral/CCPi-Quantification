cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
PROJECT (CCPiImageJPlugin)

FIND_PACKAGE(JNI REQUIRED)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

INCLUDE_DIRECTORIES(../Core)
INCLUDE_DIRECTORIES(${JNI_INCLUDE_DIRS})

ADD_LIBRARY(CCPiImageJPlugin_Cpp SHARED CCPiAccessibleVolume_wrap.cxx ../Core/CCPiAccessibleVolumeInputImages.cpp
							 ../Core/CCPiAccessibleVolumeITKImpl.cpp
							 ../Core/CCPiConsoleUserInterface.cpp)

FIND_PACKAGE(Java REQUIRED)
INCLUDE(UseJava)

FIND_JAR(ImageJ ii.jar)
set(CMAKE_JAVA_INCLUDE_PATH ${ImageJ})

ADD_JAR(CCPiImageJPlugin_Java 	AccessibleVolumePluginFilter_.java 
							CCPiAccessibleVolume.java
							CCPiAccessibleVolumeInputImages.java
							CCPiAccessibleVolumeITKImpl.java
							CCPiAccessibleVolumeJNI.java
							CCPiConsoleUserInterface.java
							CCPiImageJUserInterface.java
							CCPiUserApplicationInterface.java
							MapType.java
		OUTPUT_NAME		CCPiImageJ_Plugin.jar) 						 
							 
TARGET_LINK_LIBRARIES(CCPiImageJPlugin_Cpp ITKCommon 
							itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys)							
