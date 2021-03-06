cmake_minimum_required(VERSION 3.18)
project(fim)
set(CMAKE_CXX_STANDARD 17)

find_package(CUDA REQUIRED)
find_package(Python3 COMPONENTS Development)

include(FetchContent)

# https://github.com/nlohmann/json
FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.7.3
  GIT_SHALLOW TRUE)
FetchContent_GetProperties(json)
if(NOT json_POPULATED)
  FetchContent_Populate(json)
endif()
add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)

# https://github.com/lava/matplotlib-cpp
FetchContent_Declare(matplotlib
  GIT_REPOSITORY https://github.com/lava/matplotlib-cpp
  GIT_SHALLOW TRUE)
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
  FetchContent_Populate(matplotlib)
endif()

# https://github.com/okdshin/PicoSHA2
FetchContent_Declare(picosha2
  GIT_REPOSITORY https://github.com/okdshin/PicoSHA2
  GIT_SHALLOW TRUE)
FetchContent_GetProperties(picosha2)
if(NOT picosha2_POPULATED)
  FetchContent_Populate(picosha2)
endif()

# https://github.com/okdshin/PicoSHA2
FetchContent_Declare(threadpool
  GIT_REPOSITORY https://github.com/vit-vit/CTPL
  GIT_SHALLOW TRUE)
FetchContent_GetProperties(threadpool)
if(NOT threadpool_POPULATED)
  FetchContent_Populate(threadpool)
endif()

# https://github.com/oxyflour/ncc/blob/master/CMakeLists.txt
function(download_extract url filename)
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/build/${filename})
  file(DOWNLOAD ${url} build/${filename}.zip
    TIMEOUT 60  # seconds
    TLS_VERIFY ON)
  file(ARCHIVE_EXTRACT
    INPUT build/${filename}.zip
    DESTINATION ${filename})
endif()
endfunction(download_extract)

download_extract(
  https://prebuilt.oss-cn-shanghai.aliyuncs.com/occt-7.4-with-tbb.zip
  occt)
download_extract(
  https://prebuilt.oss-cn-shanghai.aliyuncs.com/boost_1_76_0.zip
  boost)
download_extract(
  https://prebuilt.oss-cn-shanghai.aliyuncs.com/clipper_ver6.4.2.zip
  clipper)

# generate include/generated/*.h
# https://stackoverflow.com/questions/2937128/cmake-add-custom-command-not-being-run
file(GLOB INC_HEADERS "include/*.h")
add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/include/generated/inc_chunk_cu.h
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/tool/gen.cmake
  DEPENDS src/kernel.cu tool/run.bas ${INC_HEADERS}
)
add_custom_target(
  include_generated
  DEPENDS ${CMAKE_SOURCE_DIR}/include/generated/inc_chunk_cu.h
)

link_directories(${CMAKE_SOURCE_DIR}/build/occt/occt-7.4/lib)
include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CUDA_TOOLKIT_ROOT_DIR}/include
  ${json_SOURCE_DIR}/include
  ${threadpool_SOURCE_DIR}
  ${picosha2_SOURCE_DIR}
  ${matplotlib_SOURCE_DIR}
  ${Python3_INCLUDE_DIRS}
  build/occt/occt-7.4/include
  build/boost/boost_1_76_0
  build/clipper/cpp
)

file(GLOB SOURCE_FILES "src/*.cc" "src/**/*.cc" "build/clipper/cpp/*.cpp")
list(REMOVE_ITEM SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/main.cc)
set(LINK_LIBRARIES_LIST Python3::Python nlohmann_json::nlohmann_json TKBin TKBinL TKBinTObj TKBinXCAF TKBO TKBool TKBRep TKCAF TKCDF TKD3DHost TKDCAF TKDFBrowser TKDraw TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKIVtk TKIVtkDraw TKLCAF TKMath TKMesh TKMeshVS TKOffset TKOpenGl TKPrim TKQADraw TKRWMesh TKService TKShapeView TKShHealing TKStd TKStdL TKSTEP TKSTEP209 TKSTEPAttr TKSTEPBase TKSTL TKTInspector TKTInspectorAPI TKTObj TKTObjDRAW TKToolsDraw TKTopAlgo TKTopTest TKTreeModel TKV3d TKVCAF TKView TKViewerTest TKVInspector TKVRML TKXCAF TKXDEDRAW TKXDEIGES TKXDESTEP TKXMesh TKXml TKXmlL TKXmlTObj TKXmlXCAF TKXSBase TKXSDRAW)
# Note: add following lines to your vscode settings (https://github.com/microsoft/vscode-cmake-tools/issues/1298)
#"cmake.configureSettings": {
#   "CMAKE_BUILD_TYPE": "${buildType}"
#}
if (${CMAKE_BUILD_TYPE} MATCHES Debug)
cuda_add_library(mylib ${SOURCE_FILES} "src/stl.cu" OPTIONS -D_DEBUG -Xcompiler "/wd4819,/MDd")
else ()
cuda_add_library(mylib ${SOURCE_FILES} "src/stl.cu" OPTIONS -Xcompiler "/wd4819,/MD")
endif ()

add_dependencies(mylib include_generated)

add_executable(${PROJECT_NAME} "src/main.cc")
target_link_libraries(${PROJECT_NAME} mylib ${LINK_LIBRARIES_LIST})

file(GLOB TEST_FILES "test/*.cc")
foreach(TEST_FILE ${TEST_FILES})
  file(RELATIVE_PATH TEST_NAME ${CMAKE_CURRENT_SOURCE_DIR}/test ${TEST_FILE})
  string(REPLACE ".cc" "" TEST_NAME ${TEST_NAME})
  add_executable(${TEST_NAME} ${TEST_FILE})
  target_link_libraries(${TEST_NAME} mylib ${LINK_LIBRARIES_LIST})
  add_test(${TEST_NAME} COMMAND ${TEST_NAME})
endforeach(TEST_FILE ${TEST_FILES})

IF (WIN32)
string(REPLACE "/" "\\" PROJECT_BINARY_DIR_WIN ${PROJECT_BINARY_DIR})
string(REPLACE "/" "\\" CMAKE_SOURCE_DIR_WIN ${CMAKE_SOURCE_DIR})
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD 
  COMMAND cmd /c copy 
    "${CMAKE_SOURCE_DIR_WIN}\\build\\occt\\occt-7.4\\bin\\*"
    "${PROJECT_BINARY_DIR_WIN}\\$<CONFIGURATION>\\"
  COMMENT "Copy occ binary to output directory")
ELSE()
# TODO: build for linux
ENDIF()
