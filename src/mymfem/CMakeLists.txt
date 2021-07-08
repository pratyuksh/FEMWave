target_sources(mymfem
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mybilinearform.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mymixedbf.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/local_mesh_refinement.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/utilities.cpp
)
