project(sibson)

cmake_minimum_required(VERSION 2.8.11)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(POLICY CMP0043)
    cmake_policy( SET CMP0043 OLD )
endif()

# Add specific Find.cmake modules 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules) 

find_package(CGAL)
if ( CGAL_FOUND )
include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )
include_directories (BEFORE ../../include)
else()
message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()


set(HDRS cimg.h sibson.h console.h)

set(SRCS main.cpp)

    include(AddFileDependencies)

#add_file_dependencies( main.cpp cimg.h sibson.h console.h)

find_package (Threads)
find_package(X11 REQUIRED)
include_directories(${X11_INCLUDE_DIR})

    # The executable itself.
    add_executable(sibson ${SRCS} ${HDRS})
	#set_target_properties(sibson PROPERTIES COMPILE_FLAGS "/wd4996")


    target_link_libraries( sibson  ${CMAKE_THREAD_LIBS_INIT} ${X11_LIBRARIES} )
