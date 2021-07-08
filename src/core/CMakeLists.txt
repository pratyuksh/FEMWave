target_sources(core
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/config.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/utilities.cpp
)

target_sources(mymfem
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/base_observer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/error.cpp
)

target_sources(pardiso
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/pardiso.cpp
)
