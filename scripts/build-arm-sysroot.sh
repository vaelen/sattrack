#!/bin/bash
#
# Build ARM sysroot with static libraries for cross-compilation
# Target: armv7l Linux (BeagleBone)
#
# This script downloads and cross-compiles the required dependencies
# for building SatTrack as a statically-linked ARM executable.
#
# Prerequisites:
#   sudo apt-get install gcc-arm-linux-gnueabihf g++-arm-linux-gnueabihf
#
# Usage:
#   ./scripts/build-arm-sysroot.sh
#
# The sysroot will be created at $HOME/arm-sysroot by default.
# Override with: SYSROOT=/path/to/sysroot ./scripts/build-arm-sysroot.sh

set -e

# Configuration
SYSROOT="${SYSROOT:-${HOME}/arm-sysroot}"
BUILD_DIR="${BUILD_DIR:-${HOME}/arm-sysroot-build}"
TOOLCHAIN_PREFIX="arm-linux-gnueabihf"
JOBS="${JOBS:-$(nproc)}"

# Library versions
ZLIB_VERSION="1.3.1"
OPENSSL_VERSION="3.0.13"
CURL_VERSION="8.6.0"
CURLPP_VERSION="0.8.1"
BOOST_VERSION="1.84.0"
BOOST_VERSION_UNDERSCORE="1_84_0"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check for cross-compiler
check_toolchain() {
    if ! command -v ${TOOLCHAIN_PREFIX}-gcc &> /dev/null; then
        log_error "Cross-compiler not found: ${TOOLCHAIN_PREFIX}-gcc"
        log_info "Install with: sudo apt-get install gcc-arm-linux-gnueabihf g++-arm-linux-gnueabihf"
        exit 1
    fi
    log_info "Found cross-compiler: $(${TOOLCHAIN_PREFIX}-gcc --version | head -1)"
}

# Download and extract archive
download_extract() {
    local url=$1
    local dir=$2
    local archive_name=$(basename "$url")

    if [ -d "${dir}" ]; then
        log_info "Directory ${dir} already exists, skipping download"
        return 0
    fi

    log_info "Downloading ${url}..."
    if ! wget -q --show-progress "${url}" -O "${archive_name}"; then
        log_error "Failed to download ${url}"
        exit 1
    fi

    log_info "Extracting ${archive_name}..."
    if [[ "${archive_name}" == *.tar.gz ]] || [[ "${archive_name}" == *.tgz ]]; then
        tar xzf "${archive_name}"
    elif [[ "${archive_name}" == *.tar.bz2 ]]; then
        tar xjf "${archive_name}"
    elif [[ "${archive_name}" == *.tar.xz ]]; then
        tar xJf "${archive_name}"
    else
        log_error "Unknown archive format: ${archive_name}"
        exit 1
    fi
    rm -f "${archive_name}"
}

# Build zlib (required by curl and OpenSSL)
build_zlib() {
    if [ -f "${SYSROOT}/usr/lib/libz.a" ]; then
        log_info "zlib already built, skipping"
        return 0
    fi

    log_info "=== Building zlib ${ZLIB_VERSION} ==="
    cd "${BUILD_DIR}"
    download_extract "https://zlib.net/zlib-${ZLIB_VERSION}.tar.gz" "zlib-${ZLIB_VERSION}"
    cd "zlib-${ZLIB_VERSION}"

    CC="${TOOLCHAIN_PREFIX}-gcc" \
    AR="${TOOLCHAIN_PREFIX}-ar" \
    RANLIB="${TOOLCHAIN_PREFIX}-ranlib" \
    ./configure --prefix="${SYSROOT}/usr" --static

    make -j${JOBS}
    make install
    log_info "zlib installed to ${SYSROOT}/usr"
}

# Build OpenSSL (required by curl for HTTPS)
build_openssl() {
    if [ -f "${SYSROOT}/usr/lib/libssl.a" ]; then
        log_info "OpenSSL already built, skipping"
        return 0
    fi

    log_info "=== Building OpenSSL ${OPENSSL_VERSION} ==="
    cd "${BUILD_DIR}"
    download_extract "https://www.openssl.org/source/openssl-${OPENSSL_VERSION}.tar.gz" "openssl-${OPENSSL_VERSION}"
    cd "openssl-${OPENSSL_VERSION}"

    ./Configure linux-armv4 \
        --cross-compile-prefix="${TOOLCHAIN_PREFIX}-" \
        --prefix="${SYSROOT}/usr" \
        --openssldir="${SYSROOT}/usr/ssl" \
        --with-zlib-include="${SYSROOT}/usr/include" \
        --with-zlib-lib="${SYSROOT}/usr/lib" \
        no-shared \
        no-async \
        no-dso \
        zlib

    make -j${JOBS}
    make install_sw
    log_info "OpenSSL installed to ${SYSROOT}/usr"
}

# Build libcurl (static, with minimal features)
build_curl() {
    if [ -f "${SYSROOT}/usr/lib/libcurl.a" ]; then
        log_info "libcurl already built, skipping"
        return 0
    fi

    log_info "=== Building libcurl ${CURL_VERSION} ==="
    cd "${BUILD_DIR}"
    download_extract "https://curl.se/download/curl-${CURL_VERSION}.tar.gz" "curl-${CURL_VERSION}"
    cd "curl-${CURL_VERSION}"

    # Configure with minimal features for smaller binary
    ./configure \
        --host="${TOOLCHAIN_PREFIX}" \
        --prefix="${SYSROOT}/usr" \
        --disable-shared \
        --enable-static \
        --with-openssl="${SYSROOT}/usr" \
        --with-zlib="${SYSROOT}/usr" \
        --disable-ldap \
        --disable-ldaps \
        --disable-rtsp \
        --disable-dict \
        --disable-telnet \
        --disable-tftp \
        --disable-pop3 \
        --disable-imap \
        --disable-smb \
        --disable-smtp \
        --disable-gopher \
        --disable-mqtt \
        --disable-manual \
        --disable-libcurl-option \
        --without-librtmp \
        --without-libidn2 \
        --without-libpsl \
        --without-nghttp2 \
        --without-brotli \
        --without-zstd \
        CFLAGS="-I${SYSROOT}/usr/include" \
        LDFLAGS="-L${SYSROOT}/usr/lib"

    make -j${JOBS}
    make install
    log_info "libcurl installed to ${SYSROOT}/usr"
}

# Build curlpp (C++ wrapper for curl)
build_curlpp() {
    if [ -f "${SYSROOT}/usr/lib/libcurlpp.a" ]; then
        log_info "curlpp already built, skipping"
        return 0
    fi

    log_info "=== Building curlpp ${CURLPP_VERSION} ==="
    cd "${BUILD_DIR}"

    if [ ! -d "curlpp-${CURLPP_VERSION}" ]; then
        log_info "Cloning curlpp repository..."
        git clone --depth 1 --branch v${CURLPP_VERSION} \
            https://github.com/jpbarrette/curlpp.git "curlpp-${CURLPP_VERSION}"
    fi

    cd "curlpp-${CURLPP_VERSION}"
    rm -rf build && mkdir -p build && cd build

    # Get the script directory to find the toolchain file
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    TOOLCHAIN_FILE="${SCRIPT_DIR}/../cmake/toolchains/arm-linux-gnueabihf.cmake"

    if [ ! -f "${TOOLCHAIN_FILE}" ]; then
        # Fallback: create a minimal toolchain inline
        log_warn "Toolchain file not found, using inline configuration"
        cmake .. \
            -DCMAKE_SYSTEM_NAME=Linux \
            -DCMAKE_SYSTEM_PROCESSOR=armv7l \
            -DCMAKE_C_COMPILER=${TOOLCHAIN_PREFIX}-gcc \
            -DCMAKE_CXX_COMPILER=${TOOLCHAIN_PREFIX}-g++ \
            -DCMAKE_FIND_ROOT_PATH="${SYSROOT}" \
            -DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=ONLY \
            -DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=ONLY \
            -DCMAKE_INSTALL_PREFIX="${SYSROOT}/usr" \
            -DCMAKE_PREFIX_PATH="${SYSROOT}/usr" \
            -DBUILD_SHARED_LIBS=OFF \
            -DCURL_INCLUDE_DIR="${SYSROOT}/usr/include" \
            -DCURL_LIBRARY="${SYSROOT}/usr/lib/libcurl.a"
    else
        cmake .. \
            -DCMAKE_TOOLCHAIN_FILE="${TOOLCHAIN_FILE}" \
            -DSATTRACK_ARM_SYSROOT="${SYSROOT}" \
            -DCMAKE_INSTALL_PREFIX="${SYSROOT}/usr" \
            -DCMAKE_PREFIX_PATH="${SYSROOT}/usr" \
            -DBUILD_SHARED_LIBS=OFF \
            -DCURL_INCLUDE_DIR="${SYSROOT}/usr/include" \
            -DCURL_LIBRARY="${SYSROOT}/usr/lib/libcurl.a"
    fi

    make -j${JOBS}
    make install
    log_info "curlpp installed to ${SYSROOT}/usr"
}

# Build Boost (system library only, static)
build_boost() {
    if [ -f "${SYSROOT}/usr/lib/libboost_system.a" ]; then
        log_info "Boost already built, skipping"
        return 0
    fi

    log_info "=== Building Boost ${BOOST_VERSION} ==="
    cd "${BUILD_DIR}"
    download_extract "https://boostorg.jfrog.io/artifactory/main/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION_UNDERSCORE}.tar.gz" "boost_${BOOST_VERSION_UNDERSCORE}"
    cd "boost_${BOOST_VERSION_UNDERSCORE}"

    # Bootstrap b2
    ./bootstrap.sh --prefix="${SYSROOT}/usr" --with-libraries=system

    # Create user-config.jam for cross-compilation
    cat > user-config.jam << EOF
using gcc : arm : ${TOOLCHAIN_PREFIX}-g++ ;
EOF

    # Build and install
    ./b2 install \
        --user-config=user-config.jam \
        --prefix="${SYSROOT}/usr" \
        toolset=gcc-arm \
        target-os=linux \
        architecture=arm \
        link=static \
        runtime-link=shared \
        threading=multi \
        variant=release \
        -j${JOBS}

    log_info "Boost installed to ${SYSROOT}/usr"
}

# Main
main() {
    log_info "Building ARM sysroot for SatTrack"
    log_info "Sysroot: ${SYSROOT}"
    log_info "Build dir: ${BUILD_DIR}"
    log_info "Jobs: ${JOBS}"
    echo

    check_toolchain

    # Create directories
    mkdir -p "${SYSROOT}/usr/lib" "${SYSROOT}/usr/include"
    mkdir -p "${BUILD_DIR}"

    # Build dependencies in order
    build_zlib
    build_openssl
    build_curl
    build_curlpp
    build_boost

    echo
    log_info "=== ARM sysroot build complete ==="
    log_info "Sysroot location: ${SYSROOT}"
    echo
    log_info "Static libraries built:"
    ls -la "${SYSROOT}/usr/lib/"*.a 2>/dev/null || log_warn "No .a files found"
    echo
    log_info "To build SatTrack for ARM, run:"
    log_info "  mkdir build-arm && cd build-arm"
    log_info "  cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/arm-linux-gnueabihf.cmake -DCMAKE_BUILD_TYPE=Release"
    log_info "  make -j\$(nproc)"
}

main "$@"
