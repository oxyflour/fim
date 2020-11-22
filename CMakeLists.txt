cmake_minimum_required(VERSION 3.18)
project(fim)
find_package(CUDA REQUIRED)

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

# occ

include_directories(${CMAKE_JS_INC} "${CMAKE_SOURCE_DIR}/include")
file(GLOB SOURCE_FILES "src/*.cc")
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json)
