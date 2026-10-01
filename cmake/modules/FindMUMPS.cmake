#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Optional double-precision, MPI MUMPS interface for spectral-localizer inertia.
# MUMPS_LIBRARIES may include the ordering libraries required by a static build.
find_path(
  MUMPS_INCLUDE_DIR dmumps_struc.h
  HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
  PATH_SUFFIXES include)
set(_MUMPS_REQUIRED_VARS MUMPS_INCLUDE_DIR MUMPS_LIBRARIES)
if(NOT MUMPS_LIBRARIES)
  find_library(
    MUMPS_DMUMPS_LIBRARY dmumps
    HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
    PATH_SUFFIXES lib lib64)
  find_library(
    MUMPS_COMMON_LIBRARY mumps_common
    HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
    PATH_SUFFIXES lib lib64)
  set(MUMPS_LIBRARIES ${MUMPS_DMUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY})
  find_library(
    MUMPS_PORD_LIBRARY pord
    HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
    PATH_SUFFIXES lib lib64)
  if(MUMPS_PORD_LIBRARY)
    list(APPEND MUMPS_LIBRARIES ${MUMPS_PORD_LIBRARY})
  endif()
  list(APPEND _MUMPS_REQUIRED_VARS MUMPS_DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY)
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS REQUIRED_VARS ${_MUMPS_REQUIRED_VARS})
unset(_MUMPS_REQUIRED_VARS)
if(MUMPS_FOUND AND NOT TARGET MUMPS::MUMPS)
  add_library(MUMPS::MUMPS INTERFACE IMPORTED)
  set_target_properties(
    MUMPS::MUMPS PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}"
                            INTERFACE_LINK_LIBRARIES "${MUMPS_LIBRARIES}")
endif()
mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY
                 MUMPS_PORD_LIBRARY)
