file(GLOB sources *.cpp *.c)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources new_main.cpp)

add_library(karp_core ${sources} ${headers})
target_include_directories(karp_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR})

add_executable(karp new_main.cpp)

find_package( Threads REQUIRED )
target_link_libraries(karp karp_core pthread)
add_compile_options(-W1, --no-as-needed)

if(LINK MATCHES static)
    set(BUILD_SHARED_LIBS OFF)
    set(HDF5_USE_STATIC_LIBRARIES 1)

    if (UNIX AND NOT APPLE)
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libstdc++")
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
    endif(UNIX AND NOT APPLE)

    SET_TARGET_PROPERTIES(karp karp_core PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(LINK MATCHES static)

find_package( HDF5 REQUIRED )
find_package( ZLIB REQUIRED )

if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )

if(HDF5_FOUND)
    include_directories( ${HDF5_INCLUDE_DIR} )
    target_link_libraries( karp_core ${HDF5_LIBRARIES} )
    target_link_libraries( karp ${HDF5_LIBRARIES} )
else()
    message(FATAL_ERROR "HDF5 not found. Required to output files")
endif()


if(LINK MATCHES static)
    if (UNIX AND NOT APPLE)
        target_link_libraries(karp librt.a)
    endif()
else()
    if (UNIX AND NOT APPLE)
        target_link_libraries(karp rt)
    endif()
endif(LINK MATCHES static)

install(TARGETS karp DESTINATION bin)
