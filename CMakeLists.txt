#Type 2: build main executable target

# set cmake minimum version ===================================================
cmake_minimum_required(VERSION 2.8)

# set project name ============================================================
project(ion3)

# set project options =========================================================
if(NOT DEFINED BUILD_TYPE)      # <- 1st option
#  set(BUILD_TYPE DEBUG)
  set(BUILD_TYPE RELEASE)
endif()

# clear main building variables ===============================================
set(SOURCES_LIST)
set(HEADERS_LIST)
set(INCLUDE_LIST)
set(LIBRARY_LIST)
set(DEBUG_POSTFIX)

# define target name ==========================================================
set(TARGET_NAME ion3)

# setup compilation output ====================================================
if(BUILD_TYPE MATCHES DEBUG)
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_LIST_DIR}/debug)
elseif(BUILD_TYPE MATCHES RELEASE)
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_LIST_DIR}/release)
endif()

add_subdirectory(src)

set(SOURCES_LIST ${SOURCES_LIST}
     ${CMAKE_CURRENT_LIST_DIR}/main.cpp)

# get all needble libraries ===================================================

# find_package(quadmath)

set(INCLUDE_LIST
)

set(LIBRARY_LIST
#    quadmath
)

# create & link target ========================================================
add_executable(${TARGET_NAME} ${SOURCES_LIST} ${HEADERS_LIST})
target_include_directories(${TARGET_NAME} PUBLIC ${INCLUDE_LIST})
target_link_libraries(${TARGET_NAME} ${LIBRARY_LIST})

# setup compilation definition ================================================
target_compile_definitions(${TARGET_NAME} PUBLIC
)

# setup compilation options ===================================================
set (COMMON_COMPILE_OPTIONS
#          -std=gnu++11
          -std=c++11
#          -fopenmp
#          -Wno-deprecated-declarations
)

if(BUILD_TYPE MATCHES DEBUG)
    target_compile_options(${TARGET_NAME} PUBLIC -g3 ${COMMON_COMPILE_OPTIONS})
elseif(BUILD_TYPE MATCHES RELEASE)
    target_compile_options(${TARGET_NAME} PUBLIC -O3 ${COMMON_COMPILE_OPTIONS})
endif()

