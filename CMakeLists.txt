cmake_minimum_required(VERSION 2.8)

project(karp)

#set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)

#add_compile_options(-std=c++11)
add_compile_options(-std=c++0x)
#add_compile_options(-ggdb)
add_compile_options(-O3)

if(CMAKE_BUILD_TYPE MATCHES Profile)
	add_compile_options(-g)
endif(CMAKE_BUILD_TYPE MATCHES Profile)

if(LINK MATCHES static)
    message("static build")
ELSE(LINK MATCHES shared)
    message("shared build")
ENDIF(LINK MATCHES static)

include_directories(tclap)
include_directories(karpeigen)

add_subdirectory(src)


