cmake_minimum_required(VERSION 3.18)
project(fim)
find_package(CUDA REQUIRED)
set(CMAKE_CXX_STANDARD 17)

# https://github.com/nlohmann/json
include(FetchContent)
FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.7.3)
FetchContent_GetProperties(json)
if(NOT json_POPULATED)
  FetchContent_Populate(json)
  add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
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

# https://gist.github.com/amir-saniyan/de99cee82fa9d8d615bb69f3f53b6004
function(embed_resource resource_file_name source_file_name variable_name)
  file(READ ${resource_file_name} hex_content HEX)
  string(REPEAT "[0-9a-f]" 32 column_pattern)
  string(REGEX REPLACE "(${column_pattern})" "\\1\n" content "${hex_content}")
  string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1, " content "${content}")
  string(REGEX REPLACE ", $" "" content "${content}")
  set(array_definition "static const char ${variable_name}[] =\n{\n${content}\n};")
  set(source "// Auto generated file.\n${array_definition}\n")
  file(APPEND "${source_file_name}" "${source}")
endfunction()
file(WRITE "include/generated/cst_run_bas.h" "")
embed_resource("tool/run.bas" "include/generated/cst_run_bas.h" "SRC_CST_RUN_BAS")

# https://stackoverflow.com/a/27206982
set(generated_inc "include/generated/inc_chunk_cu.h")
file(WRITE "${generated_inc}" "")
embed_resource("src/chunk.cu" "${generated_inc}" "SRC_CHUNK_CU")
file(GLOB inc_list include/*.h include/utils/*.h)
set(inc_file_list "")
foreach(inc ${inc_list})
  string(REGEX MATCH "include/(.+)$" filename ${inc})
  string(REGEX REPLACE "\\.| |-" "_" shortname "${filename}")
  string(REPLACE "/" "_" shortname "${shortname}")
  string(TOUPPER ${shortname} uppername)
  set(inc_file_list "${inc_file_list}\"${filename}\",INC_${uppername},\n")
  embed_resource(${inc} "${generated_inc}" "INC_${uppername}")
endforeach()
set(inc_file_list "static const char* INC_FILES[] = {\n${inc_file_list}\"\"\n};")
file(APPEND "${generated_inc}" "${inc_file_list}")

link_directories(${CMAKE_SOURCE_DIR}/build/occt/occt-7.4/lib)
include_directories(${CMAKE_JS_INC} "${CMAKE_SOURCE_DIR}/include" "${CUDA_TOOLKIT_ROOT_DIR}/include")
file(GLOB SOURCE_FILES "src/*.cc")
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PRIVATE build/occt/occt-7.4/include)
IF (WIN32)
target_link_libraries(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json TKBin TKBinL TKBinTObj TKBinXCAF TKBO TKBool TKBRep TKCAF TKCDF TKD3DHost TKDCAF TKDFBrowser TKDraw TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKIVtk TKIVtkDraw TKLCAF TKMath TKMesh TKMeshVS TKOffset TKOpenGl TKPrim TKQADraw TKRWMesh TKService TKShapeView TKShHealing TKStd TKStdL TKSTEP TKSTEP209 TKSTEPAttr TKSTEPBase TKSTL TKTInspector TKTInspectorAPI TKTObj TKTObjDRAW TKToolsDraw TKTopAlgo TKTopTest TKTreeModel TKV3d TKVCAF TKView TKViewerTest TKVInspector TKVRML TKXCAF TKXDEDRAW TKXDEIGES TKXDESTEP TKXMesh TKXml TKXmlL TKXmlTObj TKXmlXCAF TKXSBase TKXSDRAW)
string(REPLACE "/" "\\" PROJECT_BINARY_DIR_WIN ${PROJECT_BINARY_DIR})
string(REPLACE "/" "\\" CMAKE_SOURCE_DIR_WIN ${CMAKE_SOURCE_DIR})
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD 
  COMMAND cmd /c copy 
    "${CMAKE_SOURCE_DIR_WIN}\\build\\occt\\occt-7.4\\bin\\*"
    "${PROJECT_BINARY_DIR_WIN}\\$<CONFIGURATION>\\"
  COMMENT "Copy occ binary to output directory")
ELSE()
target_link_libraries(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json)
ENDIF()
