CMAKE_MINIMUM_REQUIRED( VERSION 3.15 )

message ("Recognized hostname:\t" $ENV{HOSTNAME})
message ("Recognized user:\t" $ENV{USER})

project(FEMWave)

# Manage dependencies
include ("cmake/GetFMT.cmake")
include ("cmake/GetJSON.cmake")
include("cmake/warning_flags.cmake")
include("cmake/clang_tidy.cmake")

# local system: nora
if (($ENV{HOSTNAME} MATCHES "nora") AND ($ENV{USER} MATCHES "prbansal"))
    
    message ("\n\n\t---------------------\t")
    message (    "\tBuild on local system\t")
    message (    "\t---------------------\t\n\n")
    
    option(USE_LIB_PARDISO "Use Pardiso library" ON)
    
    # Manage dependencies
    include ("cmake/GetGTest.cmake")
    
    # compile flags
    #set(CMAKE_CXX_COMPILER "gcc")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -O1 -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -fopenmp -DNDEBUG -msse2 -march=native")
    message ("Compiling with:\t" ${CMAKE_CXX_COMPILER})
    #message ("CXX_DEBUG flags debug:\t" ${CMAKE_CXX_FLAGS_DEBUG})
    #message ("CXX_RELEASE flags release:\t" ${CMAKE_CXX_FLAGS_RELEASE})
    
    if (CMAKE_BUILD_TYPE STREQUAL "Debug") # Debug mode
        message ("CXX_DEBUG flags debug:\t" ${CMAKE_CXX_FLAGS_DEBUG})
        add_compile_definitions(MYVERBOSE)
    endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

    if (CMAKE_BUILD_TYPE STREQUAL "Release") # Release mode
        message ("CXX_RELEASE flags release:\t" ${CMAKE_CXX_FLAGS_RELEASE})
    endif (CMAKE_BUILD_TYPE STREQUAL "Release")
    
    # external libraries
    set(PROJECT_LIB_DIR "$ENV{HOME}/Documents/softwares")
    
    # MFEM
    set(MFEM_DIR "${PROJECT_LIB_DIR}/mfem-4.1/build-serial")
    set(MFEM_LIB_DIR "${PROJECT_LIB_DIR}/mfem-4.1/build-serial")
    
    # Eigen
    set(EIGEN_INC_DIR "${PROJECT_LIB_DIR}/eigen-3.3.7")

    # Linear algebra
    message ("USE Pardiso library:\t" ${USE_LIB_PARDISO})
    add_library(linear_algebra INTERFACE)
    if(USE_LIB_PARDISO) # Pardiso library, else uses MKL Pardiso
        add_compile_definitions(LIB_PARDISO)
        include ("cmake/GetPardiso_local.cmake")
        target_link_libraries(linear_algebra INTERFACE ${PARDISO})
    endif(USE_LIB_PARDISO)
    unset(USE_LIB_PARDISO CACHE)
    target_link_libraries(linear_algebra INTERFACE -L${PROJECT_LIB_DIR}/intel_mkl_2020.0.166/install/mkl/lib/intel64 
    -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 
    -L${PROJECT_LIB_DIR}/intel_mkl_2020.0.166/install/mkl/../compiler/lib/intel64 -liomp5 -lpthread -lm -ldl)

elseif (($ENV{USER} MATCHES "prbansal"))
    
    MESSAGE("\n\n\t----------------------\t")
    MESSAGE(    "\tBuild on Euler cluster\t")
    MESSAGE(    "\t----------------------\t\n\n")
    
    option(USE_LIB_PARDISO "Use Pardiso library" OFF)
    
    # Manage dependencies
    include ("cmake/GetGTest_euler.cmake")
    
    # compile flags
    #set(CMAKE_CXX_COMPILER "gcc")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -O1")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG -msse2 -march=native")
    MESSAGE ("Compiling with:\t" ${CMAKE_CXX_COMPILER})
    #message ("CXX_DEBUG flags debug:\t" ${CMAKE_CXX_FLAGS_DEBUG})
    #message ("CXX_RELEASE flags release:\t" ${CMAKE_CXX_FLAGS_RELEASE})
    
    if (CMAKE_BUILD_TYPE STREQUAL "Debug") # Debug mode
        message ("CXX_DEBUG flags debug:\t" ${CMAKE_CXX_FLAGS_DEBUG})
        add_compile_definitions(MYVERBOSE)
    endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

    if (CMAKE_BUILD_TYPE STREQUAL "Release") # Release mode
        message ("CXX_RELEASE flags release:\t" ${CMAKE_CXX_FLAGS_RELEASE})
    endif (CMAKE_BUILD_TYPE STREQUAL "Release")
    
    # external libraries
    set(PROJECT_LIB_DIR "$ENV{HOME}/Documents/libraries")
    
    # MFEM
    set(MFEM_DIR "${PROJECT_LIB_DIR}/mfem-4.1-custom/build-serial")
    set(MFEM_LIB_DIR "${PROJECT_LIB_DIR}/mfem-4.1-custom/build-serial")
    
    # Eigen
    set(EIGEN_INC_DIR "/cluster/apps/eigen/3.2.1/x86_64/gcc_4.8.2/serial/include/eigen3")
    
     # Linear algebra
    message ("USE Pardiso library:\t" ${USE_LIB_PARDISO})
    add_library(linear_algebra INTERFACE)
    if(USE_LIB_PARDISO) # Pardiso library, else uses MKL Pardiso
        add_compile_definitions(LIB_PARDISO)
        include ("cmake/GetPardiso_euler.cmake")
        target_link_libraries(linear_algebra INTERFACE ${PARDISO})
    endif(USE_LIB_PARDISO)
    unset(USE_LIB_PARDISO CACHE)
    target_link_libraries(linear_algebra INTERFACE
    -L/cluster/apps/intel/parallel_studio_xe_2018_r1/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
    -L/cluster/apps/intel/parallel_studio_xe_2018_r1/mkl/../compiler/lib/intel64 -liomp5 -lpthread -lm -ldl)

endif()

message ("\nInclude directories:")
message ("MFEM_DIR: " ${MFEM_DIR})

message ("\nLink libraries:")
message ("MFEM_LIB_DIR: " ${MFEM_LIB_DIR})

include_directories(${MFEM_DIR}
                    ${EIGEN_INC_DIR})

link_directories(${MFEM_LIB_DIR})



##########
## Core ##
##########
add_library(core OBJECT)
target_compile_features(core PUBLIC cxx_std_17)
target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(core PUBLIC FMT)
target_link_libraries(core PUBLIC JSON)
target_link_libraries(core PUBLIC stdc++fs)
target_link_libraries(core PRIVATE warning_flags)

if(CLANG_TIDY_EXE)
    set_target_properties(core PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()


############
## MyMFEM ##
############

# MyMFEM; serial
add_library(mymfem OBJECT)
target_compile_features(mymfem PUBLIC cxx_std_17)
target_include_directories(mymfem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(mymfem PUBLIC core)
target_link_libraries(mymfem PUBLIC linear_algebra)
target_link_libraries(mymfem PRIVATE -fopenmp "${MFEM_LIB_DIR}/libmfem.a" -lmfem)

#############
## Pardiso ##
#############

add_library(pardiso OBJECT)
target_compile_features(pardiso PUBLIC cxx_std_17)
target_include_directories(pardiso PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(pardiso PUBLIC core)
target_link_libraries(pardiso PRIVATE linear_algebra)

############
## WaveO1 ##
############

# Wave; serial
add_library(waveO1Obj OBJECT)
target_link_libraries(waveO1Obj PRIVATE core)

add_executable(waveO1)
target_link_libraries(waveO1 PRIVATE core)
target_link_libraries(waveO1 PRIVATE waveO1Obj)
target_link_libraries(waveO1 PRIVATE mymfem)
target_link_libraries(waveO1 PRIVATE pardiso)


############
## WaveO2 ##
############

# Wave; serial
#add_library(waveO2Obj OBJECT)
#target_link_libraries(waveO2Obj PRIVATE core)

#add_executable(waveO2)
#target_link_libraries(waveO2 PRIVATE core)
#target_link_libraries(waveO2 PRIVATE waveO2Obj)
#target_link_libraries(waveO2 PRIVATE mymfem)
#target_link_libraries(waveO2 PRIVATE pardiso)


#############
## Testing ##
#############

# Testing, serial
add_executable(unit_tests)

target_link_libraries(unit_tests PRIVATE core)
target_link_libraries(unit_tests PRIVATE mymfem)
target_link_libraries(unit_tests PRIVATE pardiso)
target_link_libraries(unit_tests PRIVATE GTest)

#
# add sources
add_subdirectory(src)
add_subdirectory(tests)
