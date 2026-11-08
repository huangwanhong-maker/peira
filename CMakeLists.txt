# peira: designs and the mechanisms by which units come to be treated and observed, and what is identified through them.
#
# Every module is a header-only target peira::<module>. peira uses the targets of calx, morphe, physis and tyche: inside
# the infrastructure, the siblings under foundations/ (each added here when no larger project has added it). Built on its
# own, peira builds its checks (PEIRA_CHECKS); added to a larger project, only when asked.

cmake_minimum_required(VERSION 3.20)
project(peira VERSION 0.1.0 LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(PEIRA_TOP_LEVEL ON)
else()
  set(PEIRA_TOP_LEVEL OFF)
endif()
option(PEIRA_CHECKS "Build peira's checks" ${PEIRA_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# peira_sibling(<directory> <a target it defines>): the sibling library, added once.
function(peira_sibling directory target)
  if(NOT TARGET ${target})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../${directory}/CMakeLists.txt)
      add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../${directory} ${CMAKE_BINARY_DIR}/${directory})
    else()
      message(FATAL_ERROR "peira needs ${directory}: build it inside the infrastructure (foundations/${directory}).")
    endif()
  endif()
endfunction()

# The siblings in their order: calx and morphe, then physis, then tyche.
peira_sibling(calx calx::linalg)
peira_sibling(morphe morphe::topology)
peira_sibling(physis physis::units)
peira_sibling(tyche tyche::probability)
if(PEIRA_CHECKS)
  enable_testing()
endif()

# peira_module(<name> [LIBRARIES ...]): the target peira::<name> and, with the checks, check_peira_<name>, which knows where
# its reference values are (PEIRA_FIXTURES = <name>/tests/fixtures).
function(peira_module name)
  cmake_parse_arguments(ARG "" "" "LIBRARIES" ${ARGN})
  add_library(peira_${name} INTERFACE)
  add_library(peira::${name} ALIAS peira_${name})
  target_include_directories(peira_${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${name}/include)
  target_compile_features(peira_${name} INTERFACE cxx_std_20)
  target_link_libraries(peira_${name} INTERFACE ${ARG_LIBRARIES})
  if(PEIRA_CHECKS)
    calx_add_check(check_peira_${name} SOURCES ${name}/tests/check_${name}.cpp LIBRARIES peira::${name} calx::random)
    target_compile_definitions(check_peira_${name} PRIVATE PEIRA_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/${name}/tests/fixtures")
  endif()
endfunction()

# The modules in their order (decision 0003).
peira_module(design LIBRARIES calx::random tyche::probability)
peira_module(randomisation LIBRARIES peira::design tyche::probability tyche::testing)
peira_module(weighting LIBRARIES peira::design tyche::probability tyche::regression calx::linalg)
peira_module(surveys LIBRARIES peira::weighting tyche::probability)
peira_module(survival LIBRARIES tyche::probability calx::linalg)
peira_module(missing LIBRARIES peira::design tyche::probability)
