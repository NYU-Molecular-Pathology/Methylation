#!/bin/bash
## ---------------------------
## Script name: brew_installer.sh
## Purpose: Install methylation pipeline system requirements with curl and execute Rscripts
## Date Created: March 17, 2023
## Author: Jonathan Serrano
## Version: 1.1.0
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

# Check if Homebrew is installed ------------------------------------------------------------------
if ! command -v brew >/dev/null 2>&1; then
    echo -e "\nHomebrew not found. Installing Homebrew...\n"
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi

# Hardcoded URLs and function to curl download scripts --------------------------------------------
GITHUBMAIN="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"
XQUARTZURL="https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.5/"

message_curl() {
    linkname=$1
    filename=$2
    echo -e "Downloading file from Github:\n$HOME/$filename\n"
    curl -k -# -L "$linkname$filename" >"$HOME/$filename"
    chmod +rwx "$HOME/$filename"
}

# Install packages using Homebrew -----------------------------------------------------------------
install_pkgs() {
    for package in "$@"
    do
        if brew ls --versions $package > /dev/null; then
            echo -e "\n$package is already installed."
        else
            echo -e "\nInstalling $package..."
            brew install $package
        fi
    done
}

install_pkgs libmagic sqlite proj tcl-tk xz
install_pkgs aspell gdal autoconf automake gcc libgit2 openssl@3 zlib go pandoc git libffi
install_pkgs texinfo pango cairo open-mpi poppler-qt5 graphviz libopenmpt java11 zeromq libomp libtorch openjdk gmp mpfr pkg-config apache-arrow udunits mariadb-connector-c libtiff hdf5
install_pkgs llvm

brew install --cask basictex

# Check if XQuartz is installed -------------------------------------------------------------------
if [[ -d "/Applications/Utilities/XQuartz.app" ]]; then
    echo -e "XQuartz is already installed.\n"
else
    echo -e "XQuartz is not installed. Installing now...\n"
    message_curl ${XQUARTZURL} "XQuartz-2.8.5.pkg"
    sudo installer -pkg "$HOME/XQuartz-2.8.5.pkg" -target /
    echo "XQuartz installation completed."
fi

# Add JDK to PATH and link to R -------------------------------------------------------------------
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


add_pkg_path "/usr/local/opt/openjdk/bin"
add_pkg_path "/usr/local/opt/openjdk@11/bin"
sudo R CMD javareconf
R CMD config --all

# Add sqlite to PATH ------------------------------------------------------------------------------
add_pkg_path "/usr/local/opt/sqlite/bin"

# Ensure CommandLineTools is installed & Accept the Xcode license ---------------------------------
if ! xcode-select -p >/dev/null 2>&1; then
    echo -e "\nXcode Command Line Tools not found. Installing Xcode Command Line Tools...\n"
    xcode-select --install
    sudo xcode-select -s /Library/Developer/CommandLineTools
    sudo xcodebuild -license accept
fi

message_curl ${GITHUBMAIN} "all_installer.R"
Rscript --verbose "$HOME/all_installer.R"
