cmake_minimum_required(VERSION 3.18)
project(fim)
set(CMAKE_CXX_STANDARD 17)

find_package(CUDA REQUIRED)
find_package(Python3 COMPONENTS Development)

include(FetchContent)

# https://github.com/nlohmann/json
FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.7.3)
FetchContent_GetProperties(json)
if(NOT json_POPULATED)
  FetchContent_Populate(json)
  add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# https://github.com/lava/matplotlib-cpp
FetchContent_Declare(matplotlib
  GIT_REPOSITORY https://github.com/lava/matplotlib-cpp)
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
  FetchContent_Populate(matplotlib)
endif()

# https://github.com/oxyflour/ncc/blob/master/CMakeLists.txt
function(download_extract url filename)
if(NOT EXISTS ${filename})
  file(DOWNLOAD ${url} build/${filename}.zip
    TIMEOUT 60  # seconds
    TLS_VERIFY ON)
  file(ARCHIVE_EXTRACT
    INPUT build/${filename}.zip
    DESTINATION ${filename})
endif()
endfunction(download_extract)
download_extract(
  https://prebuilt.oss-cn-shanghai.aliyuncs.com/occt-7.4.zip
  occt)

# generate include/generated/*.h
add_custom_command(
  OUTPUT include_generated
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/tool/gen.cmake
)
add_custom_target(
  include_generated
  DEPENDS src/chunk.cu include/*.h include/utils/*.h)

link_directories(${CMAKE_SOURCE_DIR}/build/occt/occt-7.4/lib)
include_directories(${CMAKE_JS_INC} "${CMAKE_SOURCE_DIR}/include" "${CUDA_TOOLKIT_ROOT_DIR}/include")
file(GLOB SOURCE_FILES "src/*.cc")
add_executable(${PROJECT_NAME} ${SOURCE_FILES} include_generated)

# https://github.com/lava/matplotlib-cpp
target_include_directories(${PROJECT_NAME} PRIVATE ${matplotlib_SOURCE_DIR} build/occt/occt-7.4/include ${Python3_INCLUDE_DIRS})
IF (WIN32)
target_link_libraries(${PROJECT_NAME} PRIVATE Python3::Python nlohmann_json::nlohmann_json TKBin TKBinL TKBinTObj TKBinXCAF TKBO TKBool TKBRep TKCAF TKCDF TKD3DHost TKDCAF TKDFBrowser TKDraw TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKIVtk TKIVtkDraw TKLCAF TKMath TKMesh TKMeshVS TKOffset TKOpenGl TKPrim TKQADraw TKRWMesh TKService TKShapeView TKShHealing TKStd TKStdL TKSTEP TKSTEP209 TKSTEPAttr TKSTEPBase TKSTL TKTInspector TKTInspectorAPI TKTObj TKTObjDRAW TKToolsDraw TKTopAlgo TKTopTest TKTreeModel TKV3d TKVCAF TKView TKViewerTest TKVInspector TKVRML TKXCAF TKXDEDRAW TKXDEIGES TKXDESTEP TKXMesh TKXml TKXmlL TKXmlTObj TKXmlXCAF TKXSBase TKXSDRAW)
string(REPLACE "/" "\\" PROJECT_BINARY_DIR_WIN ${PROJECT_BINARY_DIR})
string(REPLACE "/" "\\" CMAKE_SOURCE_DIR_WIN ${CMAKE_SOURCE_DIR})
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD 
  COMMAND cmd /c copy 
    "${CMAKE_SOURCE_DIR_WIN}\\build\\occt\\occt-7.4\\bin\\*"
    "${PROJECT_BINARY_DIR_WIN}\\$<CONFIGURATION>\\"
  COMMENT "Copy occ binary to output directory")
ELSE()
# TODO
target_link_libraries(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json)
ENDIF()
