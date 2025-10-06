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

#SET(ALGOIM_HEADERS
#${gismo_externals}/algoim/quadrature_general.hpp
#)

include_directories(${gismo_externals}/libmeshb/sources)

add_library(${PROJECT_NAME} INTERFACE ${ALGOIM_SOURCES} gsLibMeshb.h)

aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

set_target_properties(${PROJECT_NAME}
PROPERTIES POSITION_INDEPENDENT_CODE ON
COMPILE_DEFINITIONS gismo_EXPORTS
FOLDER "G+Smo extensions")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "*.h")
