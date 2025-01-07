cmake_minimum_required(VERSION 2.6.2)
project(acrobot_ilqr_cpp)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread")
set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} -lX11 -ldl")

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(acrobot_ilqr_node  ${SRC_DIR}/acrobot_ilqr_test.cpp
                                  ${SRC_DIR}/acrobot_ilqr.cpp)

target_link_libraries(acrobot_ilqr_node
  ${PYTHON_LIBRARIES}
  cppad_lib
)

add_executable(acrobot_mujoco_node  ${SRC_DIR}/simulate/acrobot_mujoco.cc
                                    ${SRC_DIR}/simulate/glfw_adapter.cc
                                    ${SRC_DIR}/simulate/glfw_dispatch.cc
                                    ${SRC_DIR}/simulate/platform_ui_adapter.cc
                                    ${SRC_DIR}/simulate/simulate.cc
                                    ${SRC_DIR}/simulate/lodepng.cpp)
target_link_libraries(acrobot_mujoco_node
  ${CMAKE_DL_LIBS}
  mujoco
  glfw3
)
