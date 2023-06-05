cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(Heat)

find_package(Injection REQUIRED)
add_executable(heat "heat.cpp")
link_vnv_file(heat Heat cpp)

install(TARGETS heat DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_CURRENT_LIST_DIR}/vv-input.json DESTINATION ${CMAKE_INSTALL_PREFIX})






