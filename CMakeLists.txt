#
# fk_df cmake configuration file
#

cmake_minimum_required (VERSION 2.8)

# cmake policies
if(POLICY CMP0025) # 
  cmake_policy(SET CMP0025 NEW) # CMake 3.0 - set AppleClang instead of Clang
endif()
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW) # CMake 3.0 - macos RPATH
endif()

# load modules from cmake/ directory
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# welcome message; beginning of configuration
message (STATUS "Configuring FK dual fermion code")
project (FK_DF CXX)

# Print build type
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "git branch: ${GIT_BRANCH}")
message(STATUS "git commit hash: ${GIT_COMMIT_HASH}")

# if branch is not master - add _branch to executable names via FKPOSTFIX
if (NOT ${GIT_BRANCH} STREQUAL "master")# or ${GIT_BRANCH} strequal "Master")
    set (FKPOSTFIX "_${GIT_BRANCH}") 
    message(STATUS "Adding postfix ${FKPOSTFIX}")
endif()

# set c++-11 flags
option(AutoSetCXX11 "Auto set c++11 flags" ON)
if (AutoSetCXX11)
    include(SetCXX11) 
endif (AutoSetCXX11)

# RPATH settings (see http://www.cmake.org/Wiki/CMake_RPATH_handling) 
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
    "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# dependencies
find_package (Boost REQUIRED)
message(STATUS "Boost includes: " ${Boost_INCLUDE_DIRS} )
find_package (Eigen3 REQUIRED 3.1)
message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )
find_package (GFTools REQUIRED)
message(STATUS "GFTools includes: " ${GFTOOLS_INCLUDE_DIR} )
find_package (fftw REQUIRED)
message(STATUS "FFTW includes: " ${FFTW_INCLUDE_DIR} )

include_directories(
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${GFTOOLS_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIR}
    ./jsoncpp/include
    ./src
)

# build fk library
add_subdirectory(src)

# executables
option(Execs "Enable executables" ON)
if (Execs)
    set (dmft_lattices ${DMFT_LATTICES} CACHE STRING "Lattices for DMFT calculation" FORCE)
    set (df_lattices ${DF_LATTICES} CACHE STRING "Lattices for DF calculation" FORCE)
    add_subdirectory(prog)
    set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${PROJECT_SOURCE_DIR}/prog")
endif (Execs)

# testing
option(Testing "Enable testing" OFF)
if (Testing)
    message(STATUS "Enabled test building and running")
    add_subdirectory(tests)
    enable_testing()
endif (Testing)

# configure doxygen documentation
set(DOXYFILE_SOURCE_DIR "./src")
include(UseDoxygen)


