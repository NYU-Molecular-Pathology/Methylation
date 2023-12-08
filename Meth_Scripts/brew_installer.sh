#!/bin/bash
# Check if Homebrew is installed ------------------------------------------------------------------
if ! command -v brew >/dev/null 2>&1; then
    echo -e "\nHomebrew not found. Installing Homebrew...\n"
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi

# Install packages using Homebrew -----------------------------------------------------------------
echo -e "\nInstalling packages with brew: libmagic, sqlite, proj, tcl-tk...\n"
brew install libmagic sqlite proj tcl-tk
brew install llvm aspell gdal autoconf automake gcc libgit2 openssl@3 zlib go pandoc git libffi
brew install texinfo pango cairo open-mpi poppler-qt5 graphviz libopenmpt java11 zeromq libomp libtorch openjdk gmp mpfr pkg-config apache-arrow udunits mariadb-connector-c libtiff

# Add JDK to PATH and link to R -------------------------------------------------------------------
echo 'export PATH="/usr/local/opt/openjdk/bin:$PATH"' >> ~/.zshrc
echo 'export PATH="/usr/local/opt/openjdk/bin:$PATH"' >> ~/.bashrc
sudo R CMD javareconf

R CMD config --all

# Ensure CommandLineTools is installed ------------------------------------------------------------
if ! xcode-select -p >/dev/null 2>&1; then
    echo -e "\nXcode Command Line Tools not found. Installing Xcode Command Line Tools...\n"
    xcode-select --install
    sudo xcode-select -s /Library/Developer/CommandLineTools
    # Accept the Xcode license
    sudo xcodebuild -license accept
fi
