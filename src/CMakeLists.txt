add_subdirectory(core)
add_subdirectory(mymfem)
add_subdirectory(waveO1)
#add_subdirectory(waveO2)

target_sources(waveO1
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/waveO1.cpp
)
