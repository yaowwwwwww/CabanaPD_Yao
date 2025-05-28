# Install script for directory: /home/wuwen/program/CabanaPD/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/wuwen/program/CabanaPD/build/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/wuwen/program/CabanaPD/src/CabanaPD.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_BodyTerm.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Boundary.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Comm.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Constants.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Fields.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Force.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_ForceModels.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Geometry.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_HeatTransfer.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Input.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Integrate.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Output.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_OutputProfiles.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Particles.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Prenotch.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Solver.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Timer.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_Types.hpp"
    "/home/wuwen/program/CabanaPD/src/CabanaPD_config.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/force" TYPE FILE FILES
    "/home/wuwen/program/CabanaPD/src/force/CabanaPD_Contact.hpp"
    "/home/wuwen/program/CabanaPD/src/force/CabanaPD_Hertzian.hpp"
    "/home/wuwen/program/CabanaPD/src/force/CabanaPD_LPS.hpp"
    "/home/wuwen/program/CabanaPD/src/force/CabanaPD_PMB.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/force_models" TYPE FILE FILES
    "/home/wuwen/program/CabanaPD/src/force_models/CabanaPD_Contact.hpp"
    "/home/wuwen/program/CabanaPD/src/force_models/CabanaPD_Hertzian.hpp"
    "/home/wuwen/program/CabanaPD/src/force_models/CabanaPD_LPS.hpp"
    "/home/wuwen/program/CabanaPD/src/force_models/CabanaPD_PMB.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/home/wuwen/program/CabanaPD/src/CabanaPD_config.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/CabanaPD/CabanaPD_Targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/CabanaPD/CabanaPD_Targets.cmake"
         "/home/wuwen/program/CabanaPD/src/CMakeFiles/Export/4339873385dda2d511858f79a69bc8f8/CabanaPD_Targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/CabanaPD/CabanaPD_Targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/CabanaPD/CabanaPD_Targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/CabanaPD" TYPE FILE FILES "/home/wuwen/program/CabanaPD/src/CMakeFiles/Export/4339873385dda2d511858f79a69bc8f8/CabanaPD_Targets.cmake")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/wuwen/program/CabanaPD/src/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
