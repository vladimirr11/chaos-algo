cmake_minimum_required(VERSION 3.10.2)

project(GracefulTree LANGUAGES CXX)

set(_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}) 
set(CMAKE_CXX_STANDARD 20)

add_executable(
    ${PROJECT_NAME} 
        ${_BASE_DIR}/src/GenericTree.h
        ${_BASE_DIR}/src/GracefulTreeValidator.h
        ${_BASE_DIR}/src/TreeGenerator.h

        ${_BASE_DIR}/main.cpp
)

target_include_directories(
    ${PROJECT_NAME}
        PRIVATE 
            ${_BASE_DIR}
)
