#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

find_package(sirius 7.7.0 CONFIG REQUIRED)

# Only the Skala bridge needs the external-XC API before its numbered release.
if(CP2K_USE_LIBTORCH)
  include(CheckFortranSourceCompiles)
  include(CMakePushCheckState)
  cmake_push_check_state(RESET)
  get_target_property(_sirius_includes sirius::sirius
                      INTERFACE_INCLUDE_DIRECTORIES)
  set(CMAKE_REQUIRED_INCLUDES ${_sirius_includes})
  set(_try_compile_target_type "${CMAKE_TRY_COMPILE_TARGET_TYPE}")
  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
  unset(CP2K_SIRIUS_API_SUPPORTED CACHE)
  CHECK_Fortran_SOURCE_COMPILES(
    "program check_sirius_api
   use sirius, only: sirius_set_xc_callback, sirius_set_xc_point_callback, &
     sirius_set_xc_spinor_mode, sirius_set_xc_point_geometry_callbacks, &
     sirius_get_xc_geometry_derivatives, sirius_set_xc_point_hessian_callback, &
     sirius_set_xc_response_tolerance, sirius_set_xc_quadrature_direction_callback, &
     sirius_update_xc_point_set, sirius_check_external_xc_derivative, &
     sirius_check_external_xc_response, &
     sirius_finalize
   integer :: error
   call sirius_finalize(.false., .false., error)
   end program"
    CP2K_SIRIUS_API_SUPPORTED
    SRC_EXT
    F90)
  if(_try_compile_target_type)
    set(CMAKE_TRY_COMPILE_TARGET_TYPE "${_try_compile_target_type}")
  else()
    unset(CMAKE_TRY_COMPILE_TARGET_TYPE)
  endif()
  cmake_pop_check_state()
  unset(_try_compile_target_type)
  unset(_sirius_includes)
  if(NOT CP2K_SIRIUS_API_SUPPORTED)
    message(
      FATAL_ERROR
        "Skala with SIRIUS requires the complete external-XC and geometry-response APIs"
    )
  endif()
endif()

add_library(cp2k::sirius INTERFACE IMPORTED)
target_link_libraries(cp2k::sirius INTERFACE sirius::sirius)
get_target_property(_sirius_type sirius::sirius TYPE)
if(_sirius_type STREQUAL "SHARED_LIBRARY")
  find_package(sirius_cxx CONFIG REQUIRED)
  target_link_libraries(cp2k::sirius
                        INTERFACE "$<LINK_ONLY:sirius::sirius_cxx>")
endif()
unset(_sirius_type)
set(Sirius_FOUND TRUE)
