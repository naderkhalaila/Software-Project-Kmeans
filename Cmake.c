cmake_minimum_required(VERSION 3.20)
project(finalProject C)

set(CMAKE_C_STANDARD 90)

add_executable(finalProject spkmeans.c spkmeansmodule.c spkmeans.h)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(finalProject spkmeans.c ${PYTHON_LIBARIES})

