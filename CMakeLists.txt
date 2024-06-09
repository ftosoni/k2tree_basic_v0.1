cmake_minimum_required(VERSION 3.12)
project(k²-tree-udc)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O9 -g -DNDEBUG")

set(SOURCE_FILES basic.c bitrankw32int.c kTree.c )
set(EXECUTABLE_FILES build_tree use_tree test_tree rebuild_tree invrebuild_tree revtest_tree rebuildCheck_tree multiply)

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

#with rank

add_executable(build_tree buildk2tree.c ${SOURCE_FILES})
target_link_libraries(build_tree m)
target_compile_definitions(build_tree PRIVATE RANK_ENABLE=1)

add_executable(pagerank pagerank.c ${SOURCE_FILES} pagerank_utils.h)
target_link_libraries(pagerank m)
target_compile_definitions(pagerank PRIVATE RANK_ENABLE=1)

add_executable(pagerank_pthread pagerank_pthread.c ${SOURCE_FILES} pagerank_utils.h)
target_link_libraries(pagerank_pthread m)
target_compile_definitions(pagerank_pthread PRIVATE RANK_ENABLE=1)

#without rank

add_executable(build_tree_rd buildk2tree.c ${SOURCE_FILES})
target_link_libraries(build_tree_rd m)
target_compile_definitions(build_tree_rd PRIVATE RANK_ENABLE=0)

add_executable(pagerank_rd pagerank.c ${SOURCE_FILES} pagerank_utils.h)
target_link_libraries(pagerank_rd m)
target_compile_definitions(pagerank_rd PRIVATE RANK_ENABLE=0)

add_executable(pagerank_pthread_rd pagerank_pthread.c ${SOURCE_FILES} pagerank_utils.h)
target_link_libraries(pagerank_pthread_rd m)
target_compile_definitions(pagerank_pthread_rd PRIVATE RANK_ENABLE=0)