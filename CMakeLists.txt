cmake_minimum_required(VERSION 3.2)
project(rand-recom VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose build type (Debug|Release|RelWithDebInfo|MinSizeRel)" FORCE)
endif()

set(CMAKE_CXX_STANDARD 11)

add_definitions(-DRAND_RECOM_VERSION="${PROJECT_VERSION}")

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available

find_package(savvy REQUIRED)

add_executable(rand-recom main.cpp)

target_link_libraries(rand-recom savvy)

install(TARGETS rand-recom RUNTIME DESTINATION bin)

#set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
#set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
#set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
#include(CPack)
