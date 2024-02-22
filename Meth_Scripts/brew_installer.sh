#!/bin/bash
## ---------------------------
## Script name: brew_installer.sh
## Purpose: Install methylation pipeline system requirements with curl and execute Rscripts
## Date Created: Mar 17, 2023
## Date Last Modified: February 22, 2024
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
echo -e "\nInstalling packages with brew: libmagic, sqlite, proj, tcl-tk...\n"
brew install libmagic sqlite proj tcl-tk
brew install llvm aspell gdal autoconf automake gcc libgit2 openssl@3 zlib go pandoc git libffi
brew install texinfo pango cairo open-mpi poppler-qt5 graphviz libopenmpt java11 zeromq libomp libtorch openjdk gmp mpfr pkg-config apache-arrow udunits mariadb-connector-c libtiff hdf5
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
echo 'export PATH="/usr/local/opt/openjdk/bin:$PATH"' >> ~/.zshrc
echo 'export PATH="/usr/local/opt/openjdk/bin:$PATH"' >> ~/.bashrc
echo 'export PATH="/usr/local/opt/openjdk@11/bin:$PATH"' >> ~/.zshrc
echo 'export PATH="/usr/local/opt/openjdk@11/bin:$PATH"' >> ~/.bashrc
sudo R CMD javareconf
R CMD config --all

# Add sqlite to PATH ------------------------------------------------------------------------------
echo 'export PATH="/usr/local/opt/sqlite/bin:$PATH"' >> ~/.zshrc
echo 'export PATH="/usr/local/opt/sqlite/bin:$PATH"' >> ~/.bashrc

# Ensure CommandLineTools is installed ------------------------------------------------------------
if ! xcode-select -p >/dev/null 2>&1; then
    echo -e "\nXcode Command Line Tools not found. Installing Xcode Command Line Tools...\n"
    xcode-select --install
    sudo xcode-select -s /Library/Developer/CommandLineTools
    # Accept the Xcode license
    sudo xcodebuild -license accept
fi

message_curl ${GITHUBMAIN} "all_installer.R"
Rscript --verbose "$HOME/all_installer.R"
