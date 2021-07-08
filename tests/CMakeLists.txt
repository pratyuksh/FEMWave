target_sources(unit_tests
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_pardiso.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_mesh_point_locator.cpp
)
