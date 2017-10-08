cmake_minimum_required(VERSION 3.9)
project(c_dea)

set(CMAKE_C_STANDARD 11)

set(SOURCE_FILES main.c util.c util.h def.h lp_glpk.c lp_glpk.h)
add_executable(c_dea ${SOURCE_FILES})
target_link_libraries( c_dea glpk )