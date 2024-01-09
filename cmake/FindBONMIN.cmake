set(BONMIN_ROOT_DIR "$ENV{BONMIN_ROOT_DIR}" CACHE PATH "BONMIN root directory.")

message("Looking for Bonmin in ${BONMIN_ROOT_DIR}")

find_path(BONMIN_INCLUDE_DIR
        NAMES BonTMINLP.hpp
        HINTS /usr/local/include/coin
        HINTS /usr/local/include/coin-or
        HINTS ${BONMIN_ROOT_DIR}/include/coin
        HINTS ${BONMIN_ROOT_DIR}/include/coin-or
        HINTS ${BONMIN_ROOT_DIR}/include
)

if(APPLE)
    find_library(BONMIN_LIBRARY
            libbonmin.dylib
            HINTS /usr/local/lib
            HINTS /usr/lib
            HINTS ${BONMIN_ROOT_DIR}/lib
    )
elseif(UNIX)
    find_library(BONMIN_LIBRARY
            libbonmin.so
            HINTS /usr/local/lib
            HINTS /usr/lib
            HINTS ${BONMIN_ROOT_DIR}/lib
    )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BONMIN DEFAULT_MSG BONMIN_LIBRARY BONMIN_INCLUDE_DIR)

if(BONMIN_FOUND)
    message("—- Found Bonmin under ${BONMIN_INCLUDE_DIR}")
    set(BONMIN_INCLUDE_DIRS ${BONMIN_INCLUDE_DIR})
    set(BONMIN_LIBRARIES ${BONMIN_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(BONMIN_LIBRARIES "${BONMIN_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(BONMIN_FOUND)

mark_as_advanced(BONMIN_LIBRARY BONMIN_INCLUDE_DIR)
