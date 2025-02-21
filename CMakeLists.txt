# ======================================================================
#  dunecalib main build file
# ======================================================================


cmake_minimum_required (VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.19)

find_package(cetmodules REQUIRED)
project(dunecalib LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 10.04.03d01)

# cetbuildtools contains our cmake modules

list(APPEND CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules)
list(APPEND CMAKE_MODULE_PATH $ENV{ART_DIR}/Modules)


include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_package( art REQUIRED EXPORT )
find_package( art_root_io REQUIRED )
find_package( canvas_root_io REQUIRED )
find_package( canvas REQUIRED EXPORT )
find_package( Boost REQUIRED )
find_package(ROOT REQUIRED)
find_package( larcore REQUIRED EXPORT )
find_package( larcorealg REQUIRED EXPORT )
find_package( nuevdb REQUIRED EXPORT )
find_package( libwda REQUIRED EXPORT )
find_package( Geant4 REQUIRED EXPORT )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(dunecalib)

cet_cmake_config()
