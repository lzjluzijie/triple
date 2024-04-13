cmake_minimum_required(VERSION 3.25)

project(TRIPLE)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-pthread -Wall)
add_compile_options(-O2)
#add_compile_options(-g -ggdb)
#add_compile_options(-fsanitize=address)

file(GLOB SOURCES "src/*.cpp")
add_executable(triple ${SOURCES})

#set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_DEBUG_RUNTIME OFF)

find_package(libOTe REQUIRED
        PATHS ./libOTe
        NO_DEFAULT_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH)
target_link_libraries(triple oc::libOTe)
