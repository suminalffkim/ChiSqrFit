## CMAKE configuration file for subdirectory source
## cache the header and source files in a variable
file(GLOB SOURCE_source "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB INCLUDE_source "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
## If you want to use a root dictionary for any class in this directory, create a LinkDef.h file in include/ and comment in this line
#list(REMOVE_ITEM INCLUDE_source ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)

include_directories(include)

set(LibName source)
add_library(source ${SOURCE_source})
## If you want to use a root dictionary for any class in this directory, comment in these three lines
#set(source_LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)
#set(BASE_ROOTDIC_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
#ROOT_GENERATE_DICTIONARY(RootDictsource ${BASE_ROOTDIC_HEADER} LINKDEF ${source_LINKDEF} MODULE ${LibName})

install(TARGETS source DESTINATION lib)
## If you want to use a root dictionary for any class in this directory, comment in these this line
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${LibName}_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/lib${LibName}.rootmap DESTINATION lib)
install(FILES ${INCLUDE_source} DESTINATION include)

