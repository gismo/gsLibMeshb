######################################################################
### CMakeLists.txt --- gsLibMeshb
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

project(gsLibMeshb)

include( gsConfig)
#set(CMAKE_CXX_VISIBILITY_PRESET default)

include(gsFetch)
gismo_fetch_directory(libMeshb
  GIT_REPOSITORY https://github.com/LoicMarechal/libMeshb.git
  DESTINATION    external)

#include_directories(${gismo_externals}/libmeshb/sources)

# Set include directory
set(LIBMESHB_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/libMeshb-prefix/include CACHE INTERNAL "")

set(LIBMESHB_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/libMeshb-prefix/lib/libMeshb.7.a CACHE INTERNAL "")

# Build
include(ExternalProject)
ExternalProject_Add(libMeshb
  BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/libMeshb
  SOURCE_DIR        ${gismo_externals}/libMeshb
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
  -DCMAKE_INSTALL_MESSAGE=NEVER
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  BUILD_BYPRODUCTS "${LIBMESHB_LIBRARIES}"
  DOWNLOAD_COMMAND  ""
  UPDATE_COMMAND    ""
  )

include_directories(${gismo_SOURCE_DIR}/external/LRsplines/include;${CMAKE_CURRENT_BINARY_DIR}/LRsplines-build/include)

#Compiler wrapper
add_library(${PROJECT_NAME} OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/gsLibMeshb.h ${CMAKE_CURRENT_SOURCE_DIR}/gsLibMeshb.cpp)
# Set standard properties for all G+Smo extensions
set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS 
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions" )
  add_dependencies(${PROJECT_NAME} libMeshb)
  
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

# Add to G+Smo extensions
set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "Gismo extensions to be included")
    
# Add include directory to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${LIBMESHB_INCLUDE_DIR}
    CACHE INTERNAL "Gismo include directories")

# Link G+Smo library (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${LIBMESHB_LIBRARIES}
    CACHE INTERNAL "Gismo extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "*.h")
