target_sources(waveO1Obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_cases.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_cases_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/coefficients.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/utilities.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/assembly.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/discretisation.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/observer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/dg_error.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/dgp_error.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/energy.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/solver.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/sparse_grids_handler.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/sparse_dg_error.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/sparse_dgp_error.cpp
)
