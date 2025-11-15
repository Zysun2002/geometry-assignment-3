# Locate FREEGLUT
#
# This module defines:
#
# TO USE IT MUST SET
#
# FREEGLUT_ROOT:   Must be the folder having GL/glut.h
# or  FREEGLUT_DIR
#
# FREEGLUT_FOUND:       TRUE if the library is available
# FREEGLUT_INCLUDE_DIRS:  Directory containing header files of FREEGLUT
# FREEGLUT_LIBRARIES:     Libraries to link against when using FREEGLUT
# FREEGLUT_DLLS:          DLLs required to load the library (only on Windows)
#
# Shayan Hoshyari 21/09/2018
#
#

if(WIN32)

  find_library(FREEGLUT_LIBRARY_RELEASE
    NAMES freeglut
    PATH_SUFFIXES lib
    PATHS
    ${FREEGLUT_ROOT}
    ${FREEGLUT_DIR}
    )
  
  find_library(FREEGLUT_LIBRARY_DEBUG
    NAMES freeglutd
    PATH_SUFFIXES lib
    PATHS
    ${FREEGLUT_ROOT}
    ${FREEGLUT_DIR}
    )
  
  find_file(FREEGLUT_DLL_RELEASE
    NAMES freeglut.dll
    PATH_SUFFIXES bin
    PATHS
    ${FREEGLUT_ROOT}
    ${FREEGLUT_DIR}
    )
  
  find_file(FREEGLUT_DLL_DEBUG
    NAMES freeglutd.dll 
    PATH_SUFFIXES bin
    PATHS
    ${FREEGLUT_ROOT}
    ${FREEGLUT_DIR}
    )

  find_path(FREEGLUT_INCLUDE_DIRS GL/freeglut.h
    PATH_SUFFIXES include
    PATHS
    ${FREEGLUT_ROOT}
    ${FREEGLUT_DIR}
    )
  
  # Set up the libraries variable with optimized/debug keywords
  if(FREEGLUT_LIBRARY_DEBUG AND FREEGLUT_LIBRARY_RELEASE)
    set(FREEGLUT_LIBRARIES
      optimized ${FREEGLUT_LIBRARY_RELEASE}
      debug ${FREEGLUT_LIBRARY_DEBUG}
    )
  elseif(FREEGLUT_LIBRARY_RELEASE)
    set(FREEGLUT_LIBRARIES ${FREEGLUT_LIBRARY_RELEASE})
  elseif(FREEGLUT_LIBRARY_DEBUG)
    set(FREEGLUT_LIBRARIES ${FREEGLUT_LIBRARY_DEBUG})
  endif()
  
  # Set up the DLL variable
  if(FREEGLUT_DLL_DEBUG AND FREEGLUT_DLL_RELEASE)
    set(FREEGLUT_DLLS
      optimized ${FREEGLUT_DLL_RELEASE}
      debug ${FREEGLUT_DLL_DEBUG}
    )
  elseif(FREEGLUT_DLL_RELEASE)
    set(FREEGLUT_DLLS ${FREEGLUT_DLL_RELEASE})
  elseif(FREEGLUT_DLL_DEBUG)
    set(FREEGLUT_DLLS ${FREEGLUT_DLL_DEBUG})
  endif()

elseif(APPLE)

# COPIED FROM /usr/local/Cellar/cmake/3.14.5/share/cmake/Modules/FindGLUT.cmake
# Find glut instead!!

    find_path(FREEGLUT_INCLUDE_DIRS glut.h ${OPENGL_INCLUDE_DIR})
    find_library(GLUT_glut_LIBRARY GLUT DOC "GLUT library for OSX")
    find_library(GLUT_cocoa_LIBRARY Cocoa DOC "Cocoa framework for OSX")
    set(FREEGLUT_LIBRARIES "${GLUT_glut_LIBRARY}" "${GLUT_cocoa_LIBRARY}")
  
else()

  
  find_library(FREEGLUT_LIBRARIES
    NAMES glut  
    PATH_SUFFIXES lib
    PATHS
    ${FREEGLUT_ROOT}
    ${FREEGLUT_DIR}
    )
  
  set(FREEGLUT_DLLS "")

  find_path(FREEGLUT_INCLUDE_DIRS GL/freeglut.h
    PATH_SUFFIXES include
    PATHS
    ${FREEGLUT_ROOT}
    ${FREEGLUT_DIR}
    )

endif()


find_package(PackageHandleStandardArgs REQUIRED)
find_package_handle_standard_args(FREEGLUT
  REQUIRED_VARS FREEGLUT_INCLUDE_DIRS FREEGLUT_LIBRARIES)

