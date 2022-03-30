cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Heat)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Injection REQUIRED)
add_executable(heat "heat.cpp")
link_vnv_file(heat Heat cpp)

configure_file(${CMAKE_CURRENT_LIST_DIR}/vv-input.json ${CMAKE_CURRENT_BINARY_DIR}/vv-input.json COPYONLY) 




