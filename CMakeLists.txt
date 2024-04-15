cmake_minimum_required(VERSION 3.12)
project(TreeProject)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O9 -g -DNDEBUG")

set(SOURCE_FILES basic.c bitrankw32int.c kTree.c )
set(EXECUTABLE_FILES build_tree use_tree test_tree rebuild_tree invrebuild_tree revtest_tree rebuildCheck_tree multiply)

add_executable(build_tree buildk2tree.c ${SOURCE_FILES})
target_link_libraries(build_tree m)

add_executable(use_tree use.c ${SOURCE_FILES})
target_link_libraries(use_tree m)

add_executable(test_tree testSpeed.c ${SOURCE_FILES})
target_link_libraries(test_tree m)

add_executable(revtest_tree testSpeedPred.c ${SOURCE_FILES})
target_link_libraries(revtest_tree m)

add_executable(rebuild_tree rebuild.c ${SOURCE_FILES})
target_link_libraries(rebuild_tree m)

add_executable(rebuildCheck_tree checkrebuild.c ${SOURCE_FILES})
target_link_libraries(rebuildCheck_tree m)

add_executable(invrebuild_tree reverserebuild.c ${SOURCE_FILES})
target_link_libraries(invrebuild_tree m)

add_executable(multiply multiply.c ${SOURCE_FILES})
target_link_libraries(multiply m)

add_executable(multiply_pthread multiply_pthread.c ${SOURCE_FILES})
target_link_libraries(multiply_pthread m)
