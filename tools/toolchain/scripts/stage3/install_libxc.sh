#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxc_ver="7.1.2"
libxc_sha256="3915fac94416e4c415534223ea492ad2663f928acf27e98662c861b094a6c306"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libxc" ] && rm "${BUILDDIR}/setup_libxc"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libxc" in
  __INSTALL__)
    echo "==================== Installing LIBXC ===================="
    pkg_install_dir="${INSTALLDIR}/libxc-${libxc_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libxc-${libxc_ver} is already installed, skipping it."
    else
      retrieve_package "${libxc_sha256}" "libxc-${libxc_ver}.tar.bz2"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libxc-${libxc_ver} ] && rm -rf libxc-${libxc_ver}
      tar -xjf libxc-${libxc_ver}.tar.bz2
      cd libxc-${libxc_ver}
      # Preserve spin powers, gradient limits and physical iso-orbital derivatives.
      patch -p1 < "${SCRIPT_DIR}/stage3/libxc-${libxc_ver}-derivatives.patch" \
        > libxc_derivatives.patch.log 2>&1 || tail_excerpt libxc_derivatives.patch.log

      # Regenerate only the affected r2SCAN kernels with reproducible pure-Python
      # dependencies. Keep each functional's declared differentiability limit.
      libxc_sympy="sympy-1.13.1-py3-none-any.whl"
      libxc_sympy_sha="db36cdc64bf61b9b24578b6f7bab1ecdd2452cf008f34faa33776680c26d66f8"
      libxc_mpmath="mpmath-1.3.0-py3-none-any.whl"
      libxc_mpmath_sha="a0b2b9fe80bbcd81a6647ff13108738cfb482d481d826cc0e02f5b35e5c88d2c"
      if ! checksum "${libxc_sympy_sha}" "${BUILDDIR}/${libxc_sympy}"; then
        download_pkg_from_urlpath "${libxc_sympy_sha}" "${libxc_sympy}" \
          "https://files.pythonhosted.org/packages/b2/fe/81695a1aa331a842b582453b605175f419fe8540355886031328089d840a" \
          "${BUILDDIR}/${libxc_sympy}"
      fi
      if ! checksum "${libxc_mpmath_sha}" "${BUILDDIR}/${libxc_mpmath}"; then
        download_pkg_from_urlpath "${libxc_mpmath_sha}" "${libxc_mpmath}" \
          "https://files.pythonhosted.org/packages/43/e3/7d92a15f894aa0c9c4b49b8ee9ac9850d6e63b03c9c32c0367a13ae62209" \
          "${BUILDDIR}/${libxc_mpmath}"
      fi
      mkdir codegen-deps
      python3 -m zipfile -e "${BUILDDIR}/${libxc_sympy}" codegen-deps
      python3 -m zipfile -e "${BUILDDIR}/${libxc_mpmath}" codegen-deps
      PYTHONPATH="${PWD}/codegen-deps" python3 -s scripts/sympy2c/compile.py mgga_x_r2scan --maxorder 3 \
        > libxc_codegen.log 2>&1 || tail_excerpt libxc_codegen.log
      PYTHONPATH="${PWD}/codegen-deps" python3 -s scripts/sympy2c/compile.py mgga_c_r2scan --maxorder 3 \
        >> libxc_codegen.log 2>&1 || tail_excerpt libxc_codegen.log

      if [ "${with_gcc}" != "__DONTUSE__" ] &&
        [ "${with_intel}" = "__DONTUSE__" ] && [ "${with_amd}" = "__DONTUSE__" ]; then
        LIBXC_CFLAGS="${CFLAGS} -fno-var-tracking"
      else
        LIBXC_CFLAGS=""
      fi

      # Build exactly one variant: CUDA, HIP, or CPU.
      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        echo "Installing CUDA-only libxc into ${pkg_install_dir}"
        mkdir build-cuda
        cd build-cuda
        CFLAGS="${LIBXC_CFLAGS}" cmake \
          -DCMAKE_BUILD_TYPE="Release" \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
          -DCMAKE_INSTALL_LIBDIR="lib" \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DBUILD_SHARED_LIBS=ON \
          -DBUILD_TESTING=OFF \
          -DENABLE_FORTRAN=ON \
          -DENABLE_CUDA=ON \
          -DCMAKE_CUDA_ARCHITECTURES="${ARCH_NUM}" \
          -DMAXORDER=3 \
          .. > configure.log 2>&1 || tail_excerpt configure.log
        make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
        make install > install.log 2>&1 || tail_excerpt install.log
        cd ..
      elif [ "${ENABLE_HIP}" = "__TRUE__" ]; then
        echo "Installing HIP-only libxc into ${pkg_install_dir}"
        mkdir build-hip
        cd build-hip
        CFLAGS="${LIBXC_CFLAGS}" cmake \
          -DCMAKE_BUILD_TYPE="Release" \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
          -DCMAKE_INSTALL_LIBDIR="lib" \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DBUILD_SHARED_LIBS=ON \
          -DBUILD_TESTING=OFF \
          -DENABLE_FORTRAN=ON \
          -DENABLE_HIP=ON \
          -DCMAKE_HIP_ARCHITECTURES="${ARCH_NUM}" \
          -DMAXORDER=3 \
          .. > configure.log 2>&1 || tail_excerpt configure.log
        make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
        make install > install.log 2>&1 || tail_excerpt install.log
        cd ..
      else
        echo "Installing CPU-only libxc into ${pkg_install_dir}"
        mkdir build
        cd build
        CFLAGS="${LIBXC_CFLAGS}" cmake \
          -DCMAKE_BUILD_TYPE="Release" \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
          -DCMAKE_INSTALL_LIBDIR="lib" \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DBUILD_SHARED_LIBS=OFF \
          -DBUILD_TESTING=OFF \
          -DENABLE_FORTRAN=ON \
          -DMAXORDER=3 \
          .. > configure.log 2>&1 || tail_excerpt configure.log
        make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
        make install > install.log 2>&1 || tail_excerpt install.log
        cd ..
      fi

      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})" \
        "${SCRIPT_DIR}/stage3/libxc-${libxc_ver}-derivatives.patch"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding LIBXC from system paths ===================="
    check_lib -lxcf03 "libxc"
    check_lib -lxc "libxc"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libxc*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking LIBXC to user paths ===================="
    pkg_install_dir="${with_libxc}"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "$with_libxc" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libxc"
export LIBXC_ROOT="${pkg_install_dir}"
export LIBXC_VER="${libxc_ver}"
EOF
  if [ "$with_libxc" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libxc"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libxc" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libxc"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libxc"
