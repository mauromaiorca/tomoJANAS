cmake_minimum_required(VERSION 3.10)
project(janas)

# Set default install prefix to avoid sudo (only if not overridden)
if(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Respect the user/pip-provided prefix
else()
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Installation prefix" FORCE)
endif()

# Set build type and compiler flags
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread")

# Directories
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src/src_cpp")
set(BIN_DIR "${CMAKE_BINARY_DIR}/app_bin")

# Include directories
include_directories(${SRC_DIR})

# Create necessary directories
file(MAKE_DIRECTORY ${BIN_DIR})

# Define executables
add_executable(janas_app_meanMinMax "${SRC_DIR}/janas_app_meanMinMax.cxx")
add_executable(janas_app_starProcess "${SRC_DIR}/janas_app_starProcess.cxx")

# Set output directories
set_target_properties(janas_app_meanMinMax PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
set_target_properties(janas_app_starProcess PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

# Installation targets: install binaries to ~/.local/bin
install(TARGETS janas_app_meanMinMax janas_app_starProcess
        RUNTIME DESTINATION "bin")

# Custom install script: detect shell and update PATH using pure CMake
install(CODE [[
  # Determine shell name from environment
  string(REGEX REPLACE ".*/" "" SHELL_NAME "$ENV{SHELL}")
  if(SHELL_NAME STREQUAL "zsh")
    set(RC_FILE "$ENV{HOME}/.zshrc")
  else()
    set(RC_FILE "$ENV{HOME}/.bashrc")
  endif()

  # PATH entry to add
  set(PATH_ENTRY "export PATH=$ENV{HOME}/.local/bin:\$PATH")

  # Read existing RC file content
  file(READ "${RC_FILE}" _content)

  # Check if PATH entry is already present
  string(FIND "${_content}" "${PATH_ENTRY}" _found)
  if(_found EQUAL -1)
    # Append PATH entry
    file(APPEND "${RC_FILE}" "\n${PATH_ENTRY}\n")
  endif()
]])

message("Build complete. Run 'make install' to install without sudo and update your shell configuration.")
