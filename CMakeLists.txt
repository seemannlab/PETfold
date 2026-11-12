cmake_minimum_required(VERSION 3.18)
project(petfold VERSION 3.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard (20 or 23)")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(PETFOLD_BUILD_TOOLS "Build the command-line programs" ON)
option(PETFOLD_BUILD_PYTHON "Build the Python extension module" OFF)
option(PETFOLD_SANITIZE "Build with AddressSanitizer/UBSan" OFF)

include(GNUInstallDirs)

# ---------------------------------------------------------------- ViennaRNA
# Found via pkg-config (RNAlib2.pc); set PKG_CONFIG_PATH or VIENNARNA_ROOT.
if(DEFINED ENV{VIENNARNA_ROOT})
  list(APPEND CMAKE_PREFIX_PATH $ENV{VIENNARNA_ROOT})
  set(ENV{PKG_CONFIG_PATH} "$ENV{VIENNARNA_ROOT}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(RNALIB REQUIRED IMPORTED_TARGET RNAlib2)
find_package(Threads REQUIRED)

if(PETFOLD_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g)
  add_link_options(-fsanitize=address,undefined)
endif()

# ------------------------------------------------------------------ library
set(PFOLD_SOURCES
  src/pfold/pf_alloc.cpp src/pfold/file.cpp src/pfold/llist.cpp src/pfold/col.cpp
  src/pfold/grammar.cpp src/pfold/edouble.cpp src/pfold/matrix.cpp src/pfold/align.cpp
  src/pfold/phyl.cpp src/pfold/optimize.cpp src/pfold/inout.cpp src/pfold/newcolprob.cpp
  src/pfold/search.cpp)
set(PETFOLD_SOURCES
  src/petfold/alignment.cpp src/petfold/structure.cpp src/petfold/mea.cpp
  src/petfold/thermodynamic.cpp src/petfold/evolutionary.cpp src/petfold/petfold.cpp
  src/petfold/petcofold.cpp src/petfold/cli.cpp
  src/apps/drawdot.cpp src/apps/drawphyl.cpp)

add_library(petfold_core STATIC ${PFOLD_SOURCES} ${PETFOLD_SOURCES})
# position independent code only when it goes into the Python module
# (non-PIC code is ~10% faster for the command-line programs)
set_target_properties(petfold_core PROPERTIES POSITION_INDEPENDENT_CODE ${PETFOLD_BUILD_PYTHON})
target_include_directories(petfold_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(petfold_core PUBLIC PkgConfig::RNALIB Threads::Threads)
if(NOT PETFOLD_DATA_DIR)
  set(PETFOLD_DATA_DIR "${CMAKE_INSTALL_FULL_DATADIR}/petfold")
endif()
set_source_files_properties(src/petfold/evolutionary.cpp PROPERTIES
  COMPILE_DEFINITIONS "PETFOLD_DATA_DIR=\"${PETFOLD_DATA_DIR}\"")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(petfold_core PRIVATE -Wall)
endif()

# -------------------------------------------------------------------- tools
if(PETFOLD_BUILD_TOOLS)
  foreach(tool PETfold:petfold_main PETcofold:petcofold_main drawdot:drawdot_app drawphyl:drawphyl_app)
    string(REPLACE ":" ";" parts ${tool})
    list(GET parts 0 name)
    list(GET parts 1 src)
    add_executable(${name} src/apps/${src}.cpp)
    target_link_libraries(${name} PRIVATE petfold_core)
    install(TARGETS ${name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endforeach()
  install(FILES data/scfg.rate data/article.grm DESTINATION ${CMAKE_INSTALL_DATADIR}/petfold)
  # grammar files next to the binaries in the build tree (as in bin/ of 2.2)
  file(COPY data/scfg.rate data/article.grm DESTINATION ${CMAKE_BINARY_DIR})
endif()

# ------------------------------------------------------------------- python
if(PETFOLD_BUILD_PYTHON)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)
  # NO_EXTRAS: no LTO (RNAlib2.pc links with -fno-lto, which would leave LTO objects unlinked)
  pybind11_add_module(_core MODULE NO_EXTRAS python/src/bindings.cpp)
  target_link_libraries(_core PRIVATE petfold_core)
  set_target_properties(_core PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
  # export only PyInit__core: hide the library symbols (generic Pfold names, libRNA)
  set_target_properties(petfold_core PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_options(_core PRIVATE "-Wl,--exclude-libs,ALL")
  endif()
  install(TARGETS _core LIBRARY DESTINATION petfold)
  install(FILES data/scfg.rate data/article.grm DESTINATION petfold/data)
endif()
