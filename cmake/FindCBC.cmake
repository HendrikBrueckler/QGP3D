# - Try to find CBC
# Once done this will define
#  CBC_FOUND - System has CBC
#  CBC_INCLUDE_DIRS - The CBC include directories
#  CBC_LIBRARIES - The libraries needed to use CBC

if (CBC_INCLUDE_DIR)
  # in cache already
  set(CBC_FOUND TRUE)
  set(CBC_INCLUDE_DIRS "${CBC_INCLUDE_DIR}" )
  set(CBC_LIBRARIES "${CBC_LIBRARY}" )
else (CBC_INCLUDE_DIR)


find_path(CBC_INCLUDE_DIR
          NAMES CbcConfig.h
          PATHS "$ENV{CBC_DIR}/include/coin"
                 "/usr/include/coin"
                 "/usr/include/coin-or"
                 "/usr/local/include/coin-or"
                 "C:\\libs\\cbc\\include"
          )

find_library( CBC_LIBRARY
              NAMES Cbc libCbc
              PATHS "$ENV{CBC_DIR}/lib"
                    "/usr/lib"
                    "/usr/lib/coin"
                    "/usr/local/lib/coin-or"
                    "/usr/local/lib"
                    "C:\\libs\\cbc\\lib"
              )


set(CBC_INCLUDE_DIRS "${CBC_INCLUDE_DIR}" )
set(CBC_LIBRARIES "${CBC_LIBRARY}" )


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CBC_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CBC  DEFAULT_MSG
                                  CBC_LIBRARY CBC_INCLUDE_DIR)

mark_as_advanced(CBC_INCLUDE_DIR CBC_LIBRARY)

endif(CBC_INCLUDE_DIR)
