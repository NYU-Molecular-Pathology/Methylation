#!/bin/bash
## Script name: brew_installer.sh
## Purpose: Install methylation pipeline system requirements with curl and execute Rscripts
## Date Created: March 17, 2023
## Author: Jonathan Serrano
## Version: 1.1.0
## Copyright (c) NYULH Jonathan Serrano, 2024

# Check if Homebrew is installed --------------------------------------------
if ! command -v brew >/dev/null 2>&1; then
    echo -e "\nHomebrew not found. Installing Homebrew...\n"
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    echo >> "$HOME/.zprofile"
    echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> "$HOME/.zprofile"
    eval "$(/opt/homebrew/bin/brew shellenv)"

fi

brew update

# Hardcoded URLs and function to curl download scripts ----------------------
GITHUBMAIN="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"
XQUART_VER="XQuartz-2.8.5"
XQUARTZURL="https://github.com/XQuartz/XQuartz/releases/download/${XQUART_VER}/"

message_curl() {
    linkname=$1
    filename=$2
    echo -e "Downloading file from Github:\n$HOME/$filename\n"
    curl -k -# -L "$linkname$filename" >"$HOME/$filename"
    chmod +rwx "$HOME/$filename"
}

# Install packages using Homebrew -------------------------------------------
install_pkgs() {
    for package in "$@"
    do
        if brew ls --versions "$package" > /dev/null; then
            echo -e "$package is already installed."
        else
            echo -e "\nInstalling $package..."
            brew install "$package"
        fi
    done
}

install_pkgs libmagic sqlite proj tcl-tk xz
install_pkgs aspell gdal autoconf automake gcc libgit2 openssl@3 zlib go pandoc git libffi
install_pkgs texinfo pango cairo open-mpi poppler-qt5 graphviz libopenmpt java11 zeromq libomp pytorch openjdk@11 gmp mpfr pkg-config apache-arrow udunits mariadb-connector-c libtiff hdf5
install_pkgs llvm z3 libxml2 pcre2 bzip2 htslib blosc
install_pkgs abseil harfbuzz libtool pixman aom highway libusb pmix apache-arrow-glib hwloc libuv poppler apr icu4c libvmaf apr-util imath libvorbis popt isl libx11 portaudio jansson libxau jasper libxcb protobuf aws-sdk-cpp jpeg-turbo libxdmcp pybind11 brotli jpeg-xl libxext json-c python@3.12 c-ares jsoncpp libxrender ca-certificates krb5 libyaml qhull lame little-cms2 qt@5 cfitsio libaec rapidjson cmake libarchive lz4 re2 curl libassuan lzo readline eigen libavif m4 rtmpdump emacs libb2 make shared-mime-info epsilon libde265 sleef expat libevent minizip snappy flac mpdecimal speedtest-cli fontconfig libgcrypt freetype libgeotiff mpg123 subversion freexl netcdf swig fribidi libgpg-error netpbm gd libidn2 ninja thrift libkml npth tree-sitter gdk-pixbuf libksba nspr geos liblerc nss unbound gettext numpy unixodbc gflags libmpc openblas uriparser gh libnghttp2 utf8proc giflib libogg opencoarrays webp glib x265 glog libpng openexr wget xerces-c gnupg librsvg openldap gnutls librttopo libsndfile opus gpgme libsodium p11-kit graphite2 libspatialite zstd libssh2 grpc libtasn1 gts pinentry

# Check if basictex is installed ----------------------------------------------
if ! brew list --cask basictex &>/dev/null || command -v pdflatex > /dev/null 2>&1; then
    brew install --cask basictex
fi

# Check additional requirements for Arm64 Macs --------------------------------
if [[ $(uname -m) == "arm64" ]]; then
    brew install --cask temurin
    install_pkgs abseil gd hwloc mpdecimal pcre2 aom icu4c libgpg-error libtool pinentry gdk-pixbuf imath libheif libunistring mpg123 pixman geos isl libidn2 libusb netcdf thrift gettext jasper libkml libuv netpbm pmix tree gflags jpeg-turbo libksba libvmaf nettle udunits aws-sdk-cpp giflib jpeg-xl liblerc libvorbis npth unbound brotli json-c libx11 nspr popt unixodbc bzip2 glib jsoncpp libmpc libxau nss portaudio uriparser c-ares glog krb5 libnghttp2 libxcb numpy utf8proc ca-certificates lame libogg libxdmcp protobuf webp cfitsio gnutls libarchive openexr python@3.12 x265 eigen libassuan libpng libxrender xerces-c epsilon gpgme libavif libpq libyaml qhull xorgproto expat graphite2 libb2 librsvg little-cms2 openjpeg qt@5 flac libdap librttopo openldap rapidjson fontconfig grpc libde265 libsndfile lz4 re2 freetype gts libevent libsodium lzo opus readline freexl harfbuzz libspatialite m4 p11-kit shared-mime-info zstd fribidi libgcrypt libssh2 sleef highway libgeotiff libtasn1 minizip snappy
fi

# Check if XQuartz is installed -----------------------------------------------
if [[ -d "/Applications/Utilities/XQuartz.app" ]]; then
    echo -e "XQuartz is already installed.\n"
else
    echo -e "XQuartz is not installed. Installing now...\n"
    message_curl "${XQUARTZURL}" "${XQUART_VER}.pkg"
    sudo installer -pkg "$HOME/${XQUART_VER}.pkg" -target /
    echo "XQuartz installation completed."
fi

# Update Homebrew and upgrade outdated kegs -----------------------------------
brew doctor
brew update && brew upgrade
brew cleanup

# Add JDK to PATH and link to R -----------------------------------------------
add_pkg_path() {
    pkg_path=$1
    if ! grep -Fxq "export PATH=\"$pkg_path:\$PATH\"" ~/.zshrc
    then
        echo "export PATH=\"$pkg_path:\$PATH\"" >> ~/.zshrc
    fi

    if ! grep -Fxq "export PATH=\"$pkg_path:\$PATH\"" ~/.bashrc
    then
        echo "export PATH=\"$pkg_path:\$PATH\"" >> ~/.bashrc
    fi
}

# Add sqlite & openjdk to PATH ------------------------------------------------
OPENJDK_PATH=$(brew --prefix "openjdk@11")
SQLITE_PATH=$(brew --prefix "sqlite")
add_pkg_path "${OPENJDK_PATH}"
add_pkg_path "${SQLITE_PATH}"

# Check if the Java configuration is properly set -----------------------------
Rscript -e 'if (!requireNamespace("rJava", quietly = TRUE)) quit(JAVA_STAT = 1)'
JAVA_STAT=$?

if [ $JAVA_STAT -eq 0 ]; then 
    echo "rJava is installed"
else 
    echo "Executing: sudo R CMD javareconf"
    sudo R CMD javareconf || true
fi

# Update the PKG_CONFIG_PATH environment variable -----------------------------
ARROW_PATH=$(brew --prefix apache-arrow)
LIBOMP_PATH=$(brew --prefix libomp)

export PKG_CONFIG_PATH="$ARROW_PATH/lib/pkgconfig:$LIBOMP_PATH/lib/pkgconfig:$PKG_CONFIG_PATH"

# Ensure CommandLineTools is installed & Accept the Xcode license -------------
if ! xcode-select -p >/dev/null 2>&1; then
    echo -e "\nXcode Command Line Tools not found. Installing Xcode Command Line Tools...\n"
    xcode-select --install
    sudo xcode-select -s /Library/Developer/CommandLineTools
    sudo xcodebuild -license accept
fi

message_curl ${GITHUBMAIN} "all_installer.R"
sudo Rscript --verbose "$HOME/all_installer.R"

R CMD config --all
