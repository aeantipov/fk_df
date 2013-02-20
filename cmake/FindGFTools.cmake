# - Try to find GFTools
# Once done this will define
#  GFTOOLS_FOUND - System has LibXml2
#  GFTOOLS_INCLUDE_DIRS - The LibXml2 include directories
#  GFTOOLS_LIBRARIES - The libraries needed to use LibXml2
#  GFTOOLS_DEFINITIONS - Compiler switches required for using LibXml2

find_package(PkgConfig)
pkg_check_modules(PC_GFTOOLS QUIET gftools)
set(GFTOOLS_DEFINITIONS ${PC_GFTOOLS_CFLAGS_OTHER})
message(STATUS ${PC_GFTOOLS_INCLUDE_DIRS})

find_path(GFTOOLS_INCLUDE_DIR GridObject.hpp
          HINTS ${PC_GFTOOLS_INCLUDEDIR} ${PC_GFTOOLS_INCLUDE_DIRS}
         )
         # PATH_SUFFIXES gftools )

#find_library(GFTOOLS_LIBRARY NAMES xml2 libxml2
#             HINTS ${PC_LIBXML_LIBDIR} ${PC_LIBXML_LIBRARY_DIRS} )

#set(GFTOOLS_LIBRARIES ${GFTOOLS_LIBRARY} )
set(GFTOOLS_INCLUDE_DIRS ${GFTOOLS_INCLUDE_DIR} )

#include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GFTOOLS_FOUND to TRUE
# if all listed variables are TRUE
#find_package_handle_standard_args(GFTOOLS_INCLUDE_DIR)

mark_as_advanced(GFTOOLS_INCLUDE_DIR)
