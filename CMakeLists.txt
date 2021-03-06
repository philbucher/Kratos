project (KratosMultiphysics)
cmake_minimum_required (VERSION 2.8.6)
set (CMAKE_CXX_STANDARD 11)

# Setting some policies
# No recursive dereferencing
if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif(POLICY CMP0054)

# Uses INTERFACE_LINK_LIBRARIES instead of LINK_INTERFACE_LIBRARIES
if(POLICY CMP0022)
  cmake_policy(SET CMP0022 NEW)
endif(POLICY CMP0022)

# Find package based on <packagename>_ROOT variable
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif(POLICY CMP0074)

# Set here the version number **** only update upon tagging a release!
set (KratosMultiphysics_MAJOR_VERSION 9)
set (KratosMultiphysics_MINOR_VERSION 0)
set (KratosMultiphysics_PATCH_VERSION "Dev")

# If KRATOS_SOURCE_DIR is not defined use the CMAKE_SOURCE_DIR
if(NOT DEFINED KRATOS_SOURCE_DIR)
  message("-- No KRATOS_SOURCE_DIR is defined, using: ${CMAKE_SOURCE_DIR}")
  set (KRATOS_SOURCE_DIR "${CMAKE_SOURCE_DIR}")
endif(NOT DEFINED KRATOS_SOURCE_DIR)

if(NOT DEFINED KRATOS_BINARY_DIR)
  message("-- No KRATOS_BINARY_DIR is defined, using: ${CMAKE_BINARY_DIR}")
  set(KRATOS_BINARY_DIR ${CMAKE_BINARY_DIR})
endif(NOT DEFINED KRATOS_BINARY_DIR)

if(NOT DEFINED INSTALL_RUNKRATOS)
  message("-- No INSTALL_RUNKRATOS is defined, setting to ON")
  set(INSTALL_RUNKRATOS ON)
endif(NOT DEFINED INSTALL_RUNKRATOS)

# Link on install
if(NOT DEFINED ENV{KRATOS_INSTALL_PYTHON_USING_LINKS})
  message("-- No KRATOS_INSTALL_PYTHON_USING_LINKS is defined, setting to OFF")
  SET(INSTALL_PYTHON_USING_LINKS OFF)
else(NOT DEFINED ENV{KRATOS_INSTALL_PYTHON_USING_LINKS})
  SET(INSTALL_PYTHON_USING_LINKS $ENV{KRATOS_INSTALL_PYTHON_USING_LINKS})
endif(NOT DEFINED ENV{KRATOS_INSTALL_PYTHON_USING_LINKS})

# Setting the libs folder for the shared objects built in kratos
set(CMAKE_INSTALL_RPATH "$ORIGIN/../libs")

# If no test policy enable by default
if(NOT DEFINED KRATOS_BUILD_TESTING)
message("-- No KRATOS_BUILD_TESTING is defined. Setting to ON")
SET (KRATOS_BUILD_TESTING ON)
endif(NOT DEFINED KRATOS_BUILD_TESTING)

# If not defined mute the very verbose installation messages
if(NOT DEFINED CMAKE_INSTALL_MESSAGE)
	message("-- No CMAKE_INSTALL_MESSAGE is defined. Setting to NEVER")
	SET (CMAKE_INSTALL_MESSAGE NEVER)
endif(NOT DEFINED CMAKE_INSTALL_MESSAGE)

# Set kratos specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${KRATOS_SOURCE_DIR}/cmake_modules")

# Include cmake modules
include(DownloadLib)

# Define custom compiler build types
SET( CMAKE_CONFIGURATION_TYPES Release RelWithDebInfo Debug FullDebug Custom)
SET( BASIC_DEBUG_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
#SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" )
#SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" )
SET( CMAKE_CXX_FLAGS_DEBUG "${BASIC_DEBUG_FLAGS} -DNDEBUG -DKRATOS_DEBUG" )
SET( CMAKE_CXX_FLAGS_FULLDEBUG "${BASIC_DEBUG_FLAGS} -D_DEBUG -DKRATOS_DEBUG" )
SET( CMAKE_CXX_FLAGS_CUSTOM "${CMAKE_CXX_FLAGS_CUSTOM}" )

# Define internal CMake flags needed
SET( CMAKE_C_FLAGS_FULLDEBUG "${CMAKE_C_FLAGS_DEBUG}" )
SET( CMAKE_SHARED_LINKER_FLAGS_FULLDEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}" )
SET( CMAKE_EXE_LINKER_FLAGS_FULLDEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}" )
SET( CMAKE_C_FLAGS_CUSTOM "" )
SET( CMAKE_CXX_FLAGS_CUSTOM "" )
SET( CMAKE_SHARED_LINKER_FLAGS_CUSTOM "" )
SET( CMAKE_EXE_LINKER_FLAGS_CUSTOM "" )
SET( CMAKE_MODULE_LINKER_FLAGS_FULLDEBUG "${CMAKE_MODULE_LINKER_FLAGS}")
SET( CMAKE_MODULE_LINKER_FLAGS_CUSTOM "${CMAKE_MODULE_LINKER_FLAGS}")

# If build mode is not defined, assume Release
if(NOT CMAKE_BUILD_TYPE)
  if(DEFINED ENV{KRATOS_BUILD_TYPE})
      list(FIND CMAKE_CONFIGURATION_TYPES $ENV{KRATOS_BUILD_TYPE} FOUND_MODE)
      if(FOUND_MODE EQUAL -1)
          message("Unknown CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}', using 'Release'.")
          set(CMAKE_BUILD_TYPE Release)
      else(FOUND_MODE EQUAL -1)
          set(CMAKE_BUILD_TYPE $ENV{KRATOS_BUILD_TYPE})
      endif(FOUND_MODE EQUAL -1)
  else(DEFINED ENV{KRATOS_BUILD_TYPE})
      message("--No CMAKE_BUILD_TYPE is defined, building in 'Release' mode.")
      set(CMAKE_BUILD_TYPE Release)
  endif(DEFINED ENV{KRATOS_BUILD_TYPE})
endif(NOT CMAKE_BUILD_TYPE)
set (KratosMultiphysics_BUILD_TYPE ${CMAKE_BUILD_TYPE})


if(USE_COTIRE MATCHES ON)
  include(cotire)
  set(CMAKE_SKIP_RPATH ON)

  if(${MSVC})
    set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY True)
  endif(${MSVC})
endif(USE_COTIRE MATCHES ON)

# Get subversion data. This is done automagically by the cmakes
include (GenerateExportHeader)

# Search the SHA1 Associated with the commit in the HEAD
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    OUTPUT_VARIABLE KratosMultiphysics_SHA1_NUMBER
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else(GIT_FOUND)
  message("Git was not found on your system. SHA1 number will be set to 0.")
  set (KratosMultiphysics_SHA1_NUMBER 0)
endif(GIT_FOUND)

# Configure files depending on the build type
if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set (KratosMultiphysics_BUILD_SUFFIX "")
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set (KratosMultiphysics_BUILD_SUFFIX "-${CMAKE_BUILD_TYPE}" )
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

# Set compiler flags
if(${CMAKE_COMPILER_IS_GNUCXX})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -Wall -std=c++11")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
  endif()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -Wall")
  set(CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
  message(STATUS "additional default options were set for gcc")
  message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
endif(${CMAKE_COMPILER_IS_GNUCXX})

message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  	  if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -funroll-loops -Wall -Wno-unused-local-typedef -Wno-unknown-pragmas  ")
		set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -funroll-loops -Wall -Wno-unknown-pragmas  ")
		set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
		message("additional default options were set for clang compiler")
		message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
		message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
	  endif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	  if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -funroll-loops -Wall -Wno-unused-local-typedef -Wno-unknown-pragmas  ")
		set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops -Wall -Wno-unknown-pragmas  ")
		set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
		message("additional default options were set for clang compiler")
		message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
		message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
	  endif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
	  if(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
		if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
		  # pybind requires min. version 17, but we need at least 18:
		  message( "DEPRECATED: detected compiler as Intel " ${CMAKE_CXX_COMPILER_VERSION} )
		  message( FATAL_ERROR "Please use Version 18 or greater")
		endif()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops -lpthread -wd654 -wd10010  ")
		set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops -lpthread -wd654 -wd10010  ")
		set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
		message("additional default options were set for intel compiler")
		message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
		message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
	  endif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
else(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${CMAKE_C_COMPILER} MATCHES "icc.*$")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops  ") #-lpthread -wd654 -wd10010  ")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops  ") #-lpthread -wd654 -wd10010  ")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    message("additional default options were set for intel compiler")
    message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
    message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
  endif(${CMAKE_C_COMPILER} MATCHES "icc.*$")
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Old non-compatible versions of VS
if(${MSVC80})
  message( "DEPRECATED: detected compiler as MSVC80")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC80})

if(${MSVC90})
  message( "DEPRECATED: detected compiler as MSVC90")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC90})

if(${MSVC10})
  message( "DEPRECATED: detected compiler as MSVC10")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC10})

if(${MSVC11})
  message( "DEPRECATED: detected compiler as MSVC11")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC11})

if(${MSVC12})
  message( "DEPRECATED: detected compiler as MSVC12")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC12})

# Common flags for MSVC
if(${MSVC})
  message("Detected compiler as MSVC")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W1 /bigobj /EHsc -DBOOST_ALL_NO_LIB -D_SCL_SECURE_NO_WARNINGS")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1 /bigobj /EHsc /Zc:__cplusplus -DBOOST_ALL_NO_LIB -D_SCL_SECURE_NO_WARNINGS")
  string( REPLACE "/W3" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS} )
  string( REPLACE "/W3" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} )
endif(${MSVC})

# Specific flags for different versions of MSVC
if(${MSVC14})
  # message("Adding additional flags for MSVC14")
  # Nothing needed right now
endif(${MSVC14})

# If install prefix is set, install all the files.
SET(INSTALL_TESTING_FILES ON)
SET(INSTALL_PYTHON_FILES ON)        # To be removed when all applications are ported

#include the file with custom installation properties
include(install_function)
include(KratosDependencies)

# Logger configuration
if(KRATOS_COLORED_OUTPUT MATCHES ON)
  add_definitions(-DKRATOS_COLORED_OUTPUT)
endif(KRATOS_COLORED_OUTPUT MATCHES ON)

################### PYBIND11

# Try to use python executable from env variable
if(DEFINED ENV{PYTHON_EXECUTABLE})
  set(PYTHON_EXECUTABLE $ENV{PYTHON_EXECUTABLE})
endif(DEFINED ENV{PYTHON_EXECUTABLE})

include(pybind11Tools)

# Reset pybind11 config and remove -LTO since it gives multiple problems.
set(PYBIND11_LTO_CXX_FLAGS "" CACHE INTERNAL "")
set(PYBIND11_LTO_LINKER_FLAGS "" CACHE INTERNAL "")

# check version of Python, needs to be done after including pybind
if(${PYTHON_VERSION_MAJOR} LESS 3 OR (${PYTHON_VERSION_MAJOR} EQUAL 3 AND ${PYTHON_VERSION_MINOR} LESS 5))
    message( FATAL_ERROR "Kratos only supports Python version 3.5 and above")
endif()

# Set installation directory. TODO: Delete this and use CMAKE_INSTALL_PREFIX
if(DEFINED KRATOS_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${KRATOS_INSTALL_PREFIX} )
endif(DEFINED KRATOS_INSTALL_PREFIX)

# Set installation directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # Setting the cache path prevents it to change in case someone 'make' after modifying this file and not reconfiguring
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}" CACHE PATH "Default Install path" FORCE)
    message("-- Standard install dir ${CMAKE_INSTALL_PREFIX}")
else(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message("-- User defined install dir ${CMAKE_INSTALL_PREFIX}")
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

## Find all dependencies

# Select shared memory parallelization
if(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)
	message(FATAL_ERROR "The option \"KRATOS_EXCLUDE_OPENMP\" was removed, use KRATOS_SHARED_MEMORY_PARALLELIZATION=\"None\" instead")
endif(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)

# default to OpenMP
if(NOT DEFINED KRATOS_SHARED_MEMORY_PARALLELIZATION)
  message("\"KRATOS_SHARED_MEMORY_PARALLELIZATION\" not defined, defaulting to \"OpenMP\"")
  SET(KRATOS_SHARED_MEMORY_PARALLELIZATION "OpenMP")
endif(NOT DEFINED KRATOS_SHARED_MEMORY_PARALLELIZATION)

if (KRATOS_SHARED_MEMORY_PARALLELIZATION STREQUAL "OpenMP")
  message("Using OpenMP for shared memory parallelization")
  add_definitions( -DKRATOS_SMP_OPENMP )

	find_package(OpenMP)
	if (OPENMP_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif (OPENMP_FOUND)

elseif (KRATOS_SHARED_MEMORY_PARALLELIZATION STREQUAL "C++11")
  message("Using C++11 for shared memory parallelization")
  add_definitions( -DKRATOS_SMP_CXX11 )
  message("WARNING: Using C++11 for shared memory parallelization is highly experimental and not fully supported!")

elseif (KRATOS_SHARED_MEMORY_PARALLELIZATION STREQUAL "None")
  add_definitions( -DKRATOS_SMP_NONE )
  message("Shared memory parallelization is disabled!")

else()
  message(FATAL_ERROR "Invalid option choosen for \"KRATOS_SHARED_MEMORY_PARALLELIZATION\"! Available options are: \"OpenMP\", \"C++11\", \"None\"")
endif ()

##*****************************
# Finding and including BOOST library (version should not matter anymore)
if(DEFINED ENV{BOOST_ROOT})
  set(BOOST_ROOT $ENV{BOOST_ROOT})
endif(DEFINED ENV{BOOST_ROOT})

find_package(Boost)

if(NOT Boost_FOUND)
  # Trying to find in in external libraries
	set(BOOST_ROOT "${KRATOS_SOURCE_DIR}/external_libraries/boost")
  find_package(Boost)
endif(NOT Boost_FOUND)

set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_REALPATH ON)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

message(STATUS "Boost Include: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost Linkdir: ${Boost_LIBRARY_DIRS}")

##*****************************
# Finding cuda if needed
if(${USE_CUDA} MATCHES ON)
  find_package(CUDA QUIET REQUIRED)
endif(${USE_CUDA} MATCHES ON)


# Follow the symbolic links for the python lib only if needed
if(${INSTALL_RUNKRATOS} MATCHES ON )
    foreach( library ${PYTHON_LIBRARIES})
        get_filename_component(python_lib ${library} REALPATH)
        set(PYTHON_LIBRARIES_REALPATH ${PYTHON_LIBRARIES_REALPATH} ${python_lib})
    endforeach(library)

    set(PYTHON_LIBRARIES ${PYTHON_LIBRARIES_REALPATH} CACHE STRING "exact file of the PYTHON library" FORCE)
endif(${INSTALL_RUNKRATOS} MATCHES ON )

##*****************************

if(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
  find_package(Tetgen REQUIRED)
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON )

if(${USE_MPI} MATCHES ON )
  find_package(MPI REQUIRED)
  add_definitions( -DKRATOS_USING_MPI )
  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

  # Note: that this will be included in the compilation of the kratos core, so that it is not needed to include it somewhere else
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif(${USE_MPI} MATCHES ON )

# If we do not locally compile ZLIB we will find ZLIB in the system, if we cannot find it we will enforce the local compilation
OPTION ( FORCE_LOCAL_ZLIB_COMPILATION "Force local compilation of zlib" OFF )
if( FORCE_LOCAL_ZLIB_COMPILATION MATCHES OFF )
    # Find zlib (for gidpost)
    find_package(ZLIB)
    if( ZLIB_FOUND )
        include_directories( ${ZLIB_INCLUDE_DIRS} )
    else( ZLIB_FOUND )
        set( FORCE_LOCAL_ZLIB_COMPILATION ON )
    endif( ZLIB_FOUND )
endif( FORCE_LOCAL_ZLIB_COMPILATION MATCHES OFF )

# We local compilate ZLIB
if( FORCE_LOCAL_ZLIB_COMPILATION MATCHES ON )
  message(STATUS "Preparing local ZLIB compilation.")
  include_directories( ${KRATOS_SOURCE_DIR}/external_libraries/zlib )
  # Compile our own
  add_subdirectory( ${KRATOS_SOURCE_DIR}/external_libraries/zlib )
  # We have to add this for zconf
  include_directories( ${KRATOS_BINARY_DIR}/external_libraries/zlib )
  set( ZLIB_LIBRARIES zlib )
endif( FORCE_LOCAL_ZLIB_COMPILATION MATCHES ON )

# including the AMatrix
if( AMATRIX_DIR )
    message(STATUS "AMATRIX_DIR : ${AMATRIX_DIR}")
    include_directories("${AMATRIX_DIR}/include")
    add_definitions(-DAMATRIX_NO_IO_OPERATOR)
    add_definitions(-DKRATOS_USE_AMATRIX)
else()
    message(STATUS "AMatrix not found. Please define the AMATRIX_DIR pointing to the AMatrix root directory.")
endif()

##echo user options
message( " ")
message("CMAKE_BUILD_TYPE ........ ${CMAKE_BUILD_TYPE}")
message( " ")

# Include dir for external libraries
include_directories( ${KRATOS_SOURCE_DIR}/external_libraries )

# defines needed
add_definitions( -DKRATOS_PYTHON )

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_definitions( -fPIC )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# include subdirectories
add_subdirectory(external_libraries/gidpost)
add_subdirectory(kratos)

# Configure kratos applications
message("-- Configuring applications (ENV):")
set_property(GLOBAL PROPERTY LIST_OF_APPLICATIONS_ADDED_THROUGH_DEPENDENCIES)
if(DEFINED ENV{KRATOS_APPLICATIONS})
    foreach(APPLICATION_PATH $ENV{KRATOS_APPLICATIONS})
        get_filename_component(APPLICATION_NAME ${APPLICATION_PATH} NAME)
        list(APPEND LIST_OF_APPLICATIONS_TO_BE_COMPILED ${APPLICATION_NAME})
    endforeach(APPLICATION_PATH $ENV{KRATOS_APPLICATIONS})
endif(DEFINED ENV{KRATOS_APPLICATIONS})

if(DEFINED ENV{KRATOS_APPLICATIONS})
    foreach(APPLICATION_PATH $ENV{KRATOS_APPLICATIONS})
        get_filename_component(APPLICATION_NAME ${APPLICATION_PATH} NAME)
        if(NOT TARGET Kratos${APPLICATION_NAME})
            message("Adding application '${APPLICATION_PATH}'")
            add_subdirectory(${APPLICATION_PATH} ${CMAKE_CURRENT_BINARY_DIR}/applications/${APPLICATION_NAME})
        else(NOT TARGET Kratos${APPLICATION_NAME})
            message("[Warning] Application '${APPLICATION_PATH}' was already added")
        endif(NOT TARGET Kratos${APPLICATION_NAME})
    endforeach(APPLICATION_PATH $ENV{KRATOS_APPLICATIONS})
endif(DEFINED ENV{KRATOS_APPLICATIONS})

message("")
message("***********************************************************************")
message("")
message("        LIST OF APPLICATIONS THAT ARE GOING TO BE COMPILED:")
message("")
foreach(individual_app_name ${LIST_OF_APPLICATIONS_TO_BE_COMPILED})
  message("        ${individual_app_name}")
endforeach()
message("")
message("        APPLICATIONS ADDED TO COMPILATION THROUGH DEPENDENCIES:")
message("")
get_property(local_list GLOBAL PROPERTY LIST_OF_APPLICATIONS_ADDED_THROUGH_DEPENDENCIES)
foreach(individual_app_name ${local_list})
    if(NOT ";${LIST_OF_APPLICATIONS_TO_BE_COMPILED};" MATCHES ${individual_app_name})
        message("        ${individual_app_name}")
    endif(NOT ";${LIST_OF_APPLICATIONS_TO_BE_COMPILED};" MATCHES ${individual_app_name})
endforeach()
message("")
message("***********************************************************************")
message("")

message("List of upcoming api-breaking changes:")
message("\t Making variables non-copiable: https://github.com/KratosMultiphysics/Kratos/pull/6339")
message("\t Deprecation of NOT_FLAGS : https://github.com/KratosMultiphysics/Kratos/pull/6233")
message("\t Removing long deprecated conditions : https://github.com/KratosMultiphysics/Kratos/pull/6475")
message("\t Changing 'FindNodalNeighboursProcess' constructor signature: https://github.com/KratosMultiphysics/Kratos/pull/4947")
message("\t Minor update of import mechanism for applications: https://github.com/KratosMultiphysics/Kratos/pull/6556")
message("\t Update of Element, Condition and Scheme Interface: https://github.com/KratosMultiphysics/Kratos/issues/7381")
message("\t VariableComponent and VectorComponentAdaptor are no longer used: https://github.com/KratosMultiphysics/Kratos/pull/6686")

#compiling the triangle library
add_subdirectory(${KRATOS_SOURCE_DIR}/external_libraries/triangle)

# runkratos - compiles "runkratos and packages the python lib needed"
if(${INSTALL_RUNKRATOS} MATCHES ON )
    add_subdirectory(kratos/run_kratos)
endif(${INSTALL_RUNKRATOS} MATCHES ON )

# Clean the Module and libs install directories
install(CODE "message(STATUS \"Deleting: ${CMAKE_INSTALL_PREFIX}/KratosMultiphysics\")")
install(CODE "file(REMOVE_RECURSE \"${CMAKE_INSTALL_PREFIX}/KratosMultiphysics\")")
install(CODE "message(STATUS \"Deleting: ${CMAKE_INSTALL_PREFIX}/libs\")")
install(CODE "file(REMOVE_RECURSE \"${CMAKE_INSTALL_PREFIX}/libs\")")

# Install core files for the KratosMultiphysics python module
install(FILES "${KRATOS_SOURCE_DIR}/kratos/python_interface/__init__.py" DESTINATION KratosMultiphysics )
install(FILES "${KRATOS_SOURCE_DIR}/kratos/python_interface/kratos_globals.py" DESTINATION KratosMultiphysics )

# Install the libraries in the libs folder
install(FILES ${Boost_LIBRARIES} DESTINATION libs)
install(FILES ${EXTRA_INSTALL_LIBS} DESTINATION libs)
