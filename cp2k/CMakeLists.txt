######################################################################
# Cp2kInputExtension
# installed Avogadro library.  
# See http://avogadro.cc/ for more information.
######################################################################

# Project name - plugin name or something meaningful
project(Cp2kExtension)

# Minimum version of CMake that can be used
cmake_minimum_required(VERSION 2.6)
 
set(QT_MIN_VERSION "4.5.0")
find_package(Qt4 REQUIRED)

# Set up the build environment
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${OPENBABEL2_INCLUDE_DIR}
)

link_directories(${OPENBABEL2_LIBRARY_DIRS})

set (cp2kextension_SRCS
                cp2kextension.cpp 
                cp2kinputdialog.cpp
                cp2koutputdialog.cpp
                psfpotdialog.cpp
                )
set (cp2kextension_UIS
                    cp2kinputdialog.ui
                    cp2koutputdialog.ui
                    spectraplotdialog.ui
                    psfpotdialog.ui
                     )

avogadro_plugin(cp2kextension 
                  "${cp2kextension_SRCS}" 
                  "${cp2kextension_UIS}"   )
