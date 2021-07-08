set(JSON_BuildTests OFF CACHE INTERNAL "")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third_party/json-3.7.0)

add_library(JSON INTERFACE)
target_link_libraries(JSON INTERFACE nlohmann_json::nlohmann_json)
