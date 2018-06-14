cmake_minimum_required(VERSION 2.8)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O2")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(PkgConfig)
pkg_check_modules(SNOPT REQUIRED snopt_c snopt_cpp)

find_package(Eigen3 REQUIRED)
#find_package(SNOPT REQUIRED)
find_package(Matlab REQUIRED MX_LIBRARY)
find_package(OpenMP REQUIRED)

if (OPENMP_FOUND)
  message("OPENMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Matlab_INCLUDE_DIRS})
include_directories(${SNOPT_INCLUDE_DIRS})

link_directories(${SNOPT_LIBRARY_DIRS})

include_directories(headers)

if (CMAKE_BUILD_TYPE)
  message("The build type is " ${CMAKE_BUILD_TYPE})
  if (CMAKE_BUILD_TYPE  STREQUAL Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fkeep-inline-functions")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  endif()
endif()


message(${CMAKE_CXX_FLAGS})

#add_executable(quadrotor_ddp_est quadrotor_ddp_est.cpp)
#target_link_libraries(quadrotor_ddp_est ${SNOPT_LIBRARIES} ${Matlab_LIBRARIES} -lmat -lmx)

#add_executable(quadrotor_ddp quadrotor_ddp.cpp)
#target_link_libraries(quadrotor_ddp ${SNOPT_LIBRARIES} ${Matlab_LIBRARIES} -lmat -lmx)

#add_executable(quadrotor_ilqg quadrotor_ilqg.cpp)
#target_link_libraries(quadrotor_ilqg ${SNOPT_LIBRARIES} ${Matlab_LIBRARIES} -lmat -lmx)

#add_executable(quadrotor_lqr quadrotor_lqr.cpp)
#target_link_libraries(quadrotor_lqr ${SNOPT_LIBRARIES} ${Matlab_LIBRARIES} -lmat -lmx)

#add_executable(car_ddp car_ddp.cpp)
#target_link_libraries(car_ddp ${Matlab_LIBRARIES} -lmat -lmx)

add_executable(car car.cpp)
target_link_libraries(car ${Matlab_LIBRARIES} -lmat -lmx)

add_executable(test_boxQP test_boxQP.cc)

#add_executable(test_a test.cpp)
#target_link_libraries(test_a ${SNOPT_LIBRARIES} ${Matlab_LIBRARIES})

## #wxm
#add_executable(urdf_robot headers/urdf_robot.cpp)
#target_link_libraries(urdf_robot)
