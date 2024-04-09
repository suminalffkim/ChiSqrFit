## require cmake version 3. On centos use cmake3
cmake_minimum_required(VERSION 3.0)

## Project Name (you want to give a reasonable name here)
project(Your_Project_Name VERSION 1.0.0)

## require unix build system, since I (Henri Pekeler) did not try this anywhere else
if(NOT UNIX)
    message(FATAL_ERROR "This project must be installed on an UNIX distribution!")
endif()


## if you want you can require a specific compiler version
# if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.3)
#    message(FATAL_ERROR "GCC version must be at least 9.3! On Centos7, use devtoolset-9!")
# endif()

## set executable, library and header files for make install.
## Useful, since you can easily link with this project then.

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(ARCHIVE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(BUILD_SHARED_LIBS ON)

## You probably want to use root libraries
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

## if you want you can require specific compiler standarts
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use root compiler flags (For same c++ version)
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wextra -pthread -m64 -fno-math-errno")

## Usually, you do not need to overwrite many compiler flags but some can be useful
## If you want to specify the c++-version, you might want to use the following flags instead
# set(CMAKE_CXX_FLAGS "-std=c++11 -O2 -Wall -Wextra -pthread -m64 -fno-math-errno")
## If you use more modern compilers and want higher optimization levels, you might want to use the following flags instead
# set(CMAKE_CXX_FLAGS "-std=c++17 -O3 -Wall -Wextra -pthread -m64 -fno-math-errno")


## If you want to use any other external libraries
## In that case, do not forget to include libs also in target_link_libraries
#include_directories(<Path_To_Header_File_Directory>)
#link_directories(<Path_To_.so_File_Directory>)

## Now we can add subdirectories, where we setup our code.
## Example is given for one subdir. This is what you need to change manually!
message(STATUS "Preparing code module source")
## Have a look in the CMakeLists.txt in the sub directory for more infos on that
## and use the InitCMAKESubDirectory.sh script to create a new subdirectory. 
## Do not forget to add it here too, if you want to use it
add_subdirectory(source)

## include the header files so that the main can use them
include_directories(${CMAKE_SOURCE_DIR}/source/include)

## If you want to use a root dictionary for any class, comment in this line (See also instructions in the respective CMakeLists.txt
#include (${ROOT_USE_FILE})

## Build the executable
add_executable(Program main.cpp)

## Link the wanted libraries into the executable
## Note, that the ROOT_LIBRARIES always have to be the last think in this function call!
target_link_libraries(Program source ${ROOT_LIBRARIES})

## finaly we create a make install target
install(TARGETS Program DESTINATION bin)

## add_executable(Program2 main2.cpp)
## target_link_libraries(Program2 source ${ROOT_LIBRARIES})
## install(TARGETS Program2 DESTINATION bin)


