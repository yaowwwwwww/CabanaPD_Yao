# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/wuwen/program/CabanaPD/_deps/json-src"
  "/home/wuwen/program/CabanaPD/_deps/json-build"
  "/home/wuwen/program/CabanaPD/_deps/json-subbuild/json-populate-prefix"
  "/home/wuwen/program/CabanaPD/_deps/json-subbuild/json-populate-prefix/tmp"
  "/home/wuwen/program/CabanaPD/_deps/json-subbuild/json-populate-prefix/src/json-populate-stamp"
  "/home/wuwen/program/CabanaPD/_deps/json-subbuild/json-populate-prefix/src"
  "/home/wuwen/program/CabanaPD/_deps/json-subbuild/json-populate-prefix/src/json-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/wuwen/program/CabanaPD/_deps/json-subbuild/json-populate-prefix/src/json-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/wuwen/program/CabanaPD/_deps/json-subbuild/json-populate-prefix/src/json-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
