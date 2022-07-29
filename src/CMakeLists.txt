### Find dependencies
# mc3d
if (NOT TARGET MC3D::MC3D)
    if(EXISTS "${PROJECT_SOURCE_DIR}/extern/MC3D/CMakeLists.txt")
        add_subdirectory(${PROJECT_SOURCE_DIR}/extern/MC3D extern/MC3D EXCLUDE_FROM_ALL)
    else()
        find_package(MC3D REQUIRED)
    endif()
endif()
list(APPEND QGP3D_LIB_LIST MC3D::MC3D)

#gurobi #TODO free software alternative at some point
find_package(Gurobi REQUIRED)
list(APPEND QGP3D_LIB_LIST_PRV Gurobi::GurobiC)
list(APPEND QGP3D_LIB_LIST_PRV Gurobi::GurobiCXX)

### Add all source files
list(APPEND QGP3D_SOURCE_LIST
     "MCQuantizer.cpp"
     "ConstraintExtractor.cpp"
     "ConstraintWriter.cpp"
     "Quantizer.cpp")
list(APPEND QGP3D_HEADER_LIST
     "../include/QGP3D/PathConstraint.hpp"
     "../include/QGP3D/MCQuantizer.hpp"
     "../include/QGP3D/ConstraintExtractor.hpp"
     "../include/QGP3D/ConstraintWriter.hpp"
     "../include/QGP3D/Quantizer.hpp")

### Create target
add_library(QGP3D ${QGP3D_SOURCE_LIST} ${QGP3D_HEADER_LIST})
add_library(QGP3D::QGP3D ALIAS QGP3D)

### Configure install paths
include(CMakePackageConfigHelpers)
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR "CMake")
else()
  set(DEF_INSTALL_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/QGP3D")
endif()
set(QGP3D_INSTALL_CMAKE_DIR "${DEF_INSTALL_CMAKE_DIR}" CACHE PATH "Installation directory for CMake files")
set(QGP3D_INSTALL_LIB_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "Installation directory for libraries")
set(QGP3D_INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Installation directory for header files")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../include")
set(CONF_LIB_DIRS "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
configure_package_config_file(../QGP3DConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/QGP3DConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
# ... for the install tree
set(CONF_INCLUDE_DIRS "${QGP3D_INSTALL_INCLUDE_DIR}")
set(CONF_LIB_DIRS "${QGP3D_INSTALL_LIB_DIR}")
configure_package_config_file(../QGP3DConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/QGP3DConfig.cmake"
  INSTALL_DESTINATION "${QGP3D_INSTALL_CMAKE_DIR}")

### Include own headers with public access
target_include_directories(QGP3D
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>"
                           "$<INSTALL_INTERFACE:${QGP3D_INSTALL_INCLUDE_DIR}>")

### Link with dependencies
target_link_libraries_system(QGP3D PUBLIC ${QGP3D_LIB_LIST})
target_link_libraries_system(QGP3D PRIVATE ${QGP3D_LIB_LIST_PRV})

### Set properties
# properties
set_target_properties(QGP3D PROPERTIES
                      DEBUG_POSTFIX "d"
                      VERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}
                      SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}
                      PUBLIC_HEADER "${QGP3D_HEADER_LIST}"
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED ON)
# compile options
list(APPEND QGP3D_COMPILE_OPTIONS_PRV "-Wall" "-Wextra" "-Wshadow" "-Wunused" "-Wunused-function")
target_compile_options(QGP3D PRIVATE ${QGP3D_COMPILE_OPTIONS_PRV})
# preprocessor defines
if (NOT QGP3D_ENABLE_LOGGING)
    list(APPEND QGP3D_COMPILE_DEFINITIONS_PRV "GOOGLE_STRIP_LOG=10")
endif()
target_compile_definitions(QGP3D PRIVATE ${QGP3D_COMPILE_DEFINITIONS_PRV})

### Install and export
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/QGP3DConfig.cmake"
        DESTINATION "${QGP3D_INSTALL_CMAKE_DIR}")

install(TARGETS       QGP3D
        EXPORT        QGP3DTargets
        LIBRARY       DESTINATION "${QGP3D_INSTALL_LIB_DIR}"
        ARCHIVE       DESTINATION "${QGP3D_INSTALL_LIB_DIR}"
        PUBLIC_HEADER DESTINATION "${QGP3D_INSTALL_INCLUDE_DIR}/QGP3D")

install(EXPORT QGP3DTargets
        FILE QGP3DTargets.cmake
        NAMESPACE QGP3D::
        DESTINATION "${QGP3D_INSTALL_CMAKE_DIR}")

export(EXPORT QGP3DTargets
       NAMESPACE QGP3D::)
