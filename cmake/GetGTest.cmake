find_package(Threads REQUIRED)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest-release-1.8.1)

add_library(GTest INTERFACE)

target_link_libraries(
  GTest INTERFACE
  gtest$<$<AND:$<OR:$<PLATFORM_ID:Windows>,$<AND:$<PLATFORM_ID:Darwin>,$<CXX_COMPILER_ID:AppleClang>>>,$<CONFIG:Debug>>:d>
  ${THREADS_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  )
