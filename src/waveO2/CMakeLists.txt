target_sources(waveO2Obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_cases.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_cases_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/coefficients.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/utilities.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/assembly.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/discretisation.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/observer.cpp
  #PRIVATE ${CMAKE_CURRENT_LIST_DIR}/solver.cpp
  #PRIVATE ${CMAKE_CURRENT_LIST_DIR}/sparse_grids_handler.cpp
)
