#!/bin/bash
## ---------------------------
## Script name: brew_installer.sh
## Purpose: Install methylation pipeline system requirements with curl and execute Rscripts
## Date Created: March 17, 2023
## Author: Jonathan Serrano
## Version: 1.1.0
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

# Check if Homebrew is installed --------------------------------------------
if ! command -v brew >/dev/null 2>&1; then
    echo -e "\nHomebrew not found. Installing Homebrew...\n"
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi

# Hardcoded URLs and function to curl download scripts ----------------------
GITHUBMAIN="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"
XQUARTZURL="https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.5/"

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
        if brew ls --versions $package > /dev/null; then
            echo -e "$package is already installed."
        else
            echo -e "\nInstalling $package..."
            brew install $package
        fi
    done
}

install_pkgs libmagic sqlite proj tcl-tk xz
install_pkgs aspell gdal autoconf automake gcc libgit2 openssl@3 zlib go pandoc git libffi
install_pkgs texinfo pango cairo open-mpi poppler-qt5 graphviz libopenmpt java11 zeromq libomp pytorch openjdk gmp mpfr pkg-config apache-arrow udunits mariadb-connector-c libtiff hdf5
install_pkgs llvm z3 libxml2

# Check if basictex is installed ----------------------------------------------
if ! brew list --cask basictex &>/dev/null || command -v pdflatex > /dev/null 2>&1; then
    brew install --cask basictex
fi

# Check additional requirements for Arm64 Macs --------------------------------
if [[ $(uname -m) == "arm64" ]]; then
    brew install --cask adoptopenjdk
    install_pkgs abseil gd hwloc mpdecimal pcre2 aom icu4c libgpg-error libtool pinentry apache-arrow gdk-pixbuf imath libheif libunistring mpg123 pixman geos isl libidn2 libusb netcdf thrift gettext jasper libkml libuv netpbm pmix tree gflags jpeg-turbo libksba libvmaf nettle udunits aws-sdk-cpp giflib jpeg-xl liblerc libvorbis npth unbound brotli json-c libx11 nspr popt unixodbc bzip2 glib jsoncpp libmpc libxau nss portaudio uriparser c-ares glog krb5 libnghttp2 libxcb numpy utf8proc ca-certificates lame libogg libxdmcp protobuf webp cfitsio gnutls libarchive openexr python@3.12 x265 eigen libassuan libpng libxrender xerces-c epsilon gpgme libavif libpq libyaml openjdk@11 qhull xorgproto expat graphite2 libb2 librsvg little-cms2 openjpeg qt@5 flac libdap librttopo openldap rapidjson fontconfig grpc libde265 libsndfile lz4 re2 freetype gts libevent libsodium lzo opus readline freexl harfbuzz libspatialite m4 p11-kit shared-mime-info zstd fribidi libgcrypt libssh2 sleef highway libgeotiff libtasn1 minizip snappy
fi

# Check if XQuartz is installed -----------------------------------------------
if [[ -d "/Applications/Utilities/XQuartz.app" ]]; then
    echo -e "XQuartz is already installed.\n"
else
    echo -e "XQuartz is not installed. Installing now...\n"
    message_curl ${XQUARTZURL} "XQuartz-2.8.5.pkg"
    sudo installer -pkg "$HOME/XQuartz-2.8.5.pkg" -target /
    echo "XQuartz installation completed."
fi

# Update Homebrew and upgrade outdated kegs -----------------------------------
brew doctor
brew update && brew upgrade
# brew list --formula | while read -r keg; do
#     if ! brew link --overwrite --force "$keg" 2>/dev/null; then
#         brew unlink "$keg" && brew link --overwrite --force "$keg"
#     fi
# done
brew cleanup

# Add JDK to PATH and link to R -----------------------------------------------
add_pkg_path() {
    path=$1
    for file in ~/.zshrc ~/.bashrc
    do
        if ! grep -Fxq "export PATH=\"$path:\$PATH\"" $file
        then
            echo "export PATH=\"$path:\$PATH\"" >> $file
        fi
    done
}

# Add sqlite & openjdk to PATH ------------------------------------------------
add_pkg_path "/usr/local/opt/sqlite/bin"
add_pkg_path "/usr/local/opt/openjdk/bin"

# Check if the Java configuration is properly set -----------------------------
# if ! R CMD javareconf -e | grep -q "Java interpreter : /usr/bin/java"; then
#     sudo R CMD javareconf
# fi

# Update the PKG_CONFIG_PATH environment variable -----------------------------
ARROW_PATH=$(brew --prefix apache-arrow)
OPENJDK_PATH=$(brew --prefix openjdk)
LIBOMP_PATH=$(brew --prefix libomp)

export PKG_CONFIG_PATH="$ARROW_PATH/lib/pkgconfig:$OPENJDK_PATH/lib/pkgconfig:$LIBOMP_PATH/lib/pkgconfig:$PKG_CONFIG_PATH"

# Ensure CommandLineTools is installed & Accept the Xcode license -------------
if ! xcode-select -p >/dev/null 2>&1; then
    echo -e "\nXcode Command Line Tools not found. Installing Xcode Command Line Tools...\n"
    xcode-select --install
    sudo xcode-select -s /Library/Developer/CommandLineTools
    sudo xcodebuild -license accept
fi

#R CMD config --all

message_curl ${GITHUBMAIN} "all_installer.R"
Rscript --verbose "$HOME/all_installer.R"
