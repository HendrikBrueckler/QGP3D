### Better safe than sorry...
list(APPEND QGP3D_COMPILE_OPTIONS "-Wall" "-Wextra" "-Werror" "-Wshadow")

### Add all source files
list(APPEND QGP3D_SOURCE_LIST
     "MCQuantizer.cpp"
     "ConstraintExtractor.cpp"
     "ConstraintWriter.cpp"
     "Quantizer.cpp"
     )
add_library(${QGP3D_LIB_NAME} STATIC ${QGP3D_SOURCE_LIST})

### Find dependencies
# mc3d
if (NOT TARGET mc3d)
    add_subdirectory(${PROJECT_SOURCE_DIR}/extern/MC3D extern/MC3D EXCLUDE_FROM_ALL)
endif()
list(APPEND QGP3D_LIB_LIST mc3d)

#gurobi #TODO free software alternative at some point
find_package(GUROBI REQUIRED)
list(APPEND QGP3D_LIB_LIST Gurobi::GurobiC)
list(APPEND QGP3D_LIB_LIST Gurobi::GurobiCXX)

### Link with dependencies
target_link_libraries_system(${QGP3D_LIB_NAME} PUBLIC OpenVolumeMesh::OpenVolumeMesh ${QGP3D_LIB_LIST})

### Set compile options
set_target_properties(${QGP3D_LIB_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_compile_options(${QGP3D_LIB_NAME} PRIVATE ${QGP3D_COMPILE_OPTIONS})
target_compile_definitions(${QGP3D_LIB_NAME} PRIVATE ${QGP3D_COMPILE_DEFINITIONS})

### Include own headers with public access
target_include_directories(${QGP3D_LIB_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

### Include headers as system headers to prevent warnings/errors from strict compile definitions
target_include_directories(${QGP3D_LIB_NAME} SYSTEM PUBLIC ${GMPXX_INCLUDE_DIR})
target_include_directories(${QGP3D_LIB_NAME} SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/MC3D/include)
target_include_directories(${QGP3D_LIB_NAME} SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/MC3D/extern/OpenVolumeMesh/src)
target_include_directories(${QGP3D_LIB_NAME} SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/MC3D/extern/glog/src)
