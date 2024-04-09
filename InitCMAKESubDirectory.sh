## Script to initiate a subdirectory for a cmake projekt
if [ $# -eq 0 ]
 then
  echo "Please give an argument to this script!"
  echo "(The name of the subdirectory you want to initialize)"
  exit 1
fi

mkdir $1
echo "## CMAKE configuration file for subdirectory $1
## cache the header and source files in a variable
file(GLOB SOURCE_$1 \"\${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp\")
file(GLOB INCLUDE_$1 \"\${CMAKE_CURRENT_SOURCE_DIR}/include/*.h\")
## If you want to use a root dictionary for any class in this directory, create a LinkDef.h file in include/ and comment in this line
#list(REMOVE_ITEM INCLUDE_$1 \${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)

include_directories(include)

set(LibName $1)
add_library($1 \${SOURCE_$1})
## If you want to use a root dictionary for any class in this directory, comment in these three lines
#set($1_LINKDEF \${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)
#set(BASE_ROOTDIC_HEADER \${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
#ROOT_GENERATE_DICTIONARY(RootDict$1 \${BASE_ROOTDIC_HEADER} LINKDEF \${$1_LINKDEF} MODULE \${LibName})

install(TARGETS $1 DESTINATION lib)
## If you want to use a root dictionary for any class in this directory, comment in these this line
#install(FILES \${CMAKE_CURRENT_BINARY_DIR}/lib\${LibName}_rdict.pcm \${CMAKE_CURRENT_BINARY_DIR}/lib\${LibName}.rootmap DESTINATION lib)
install(FILES \${INCLUDE_$1} DESTINATION include)
" > $1/CMakeLists.txt

mkdir $1/include
mkdir $1/src

echo "#pragma once
class A {
    public:
    A();
};" > $1/include/Example.h

echo "#include \"Example.h\"

#include <iostream>

A::A() { std::cout << \"Created Object of Class A\" << std::endl; } " > $1/src/Example.cpp

echo "Do not forget to add sub directory $1 to your CMakeLists.txt!"

