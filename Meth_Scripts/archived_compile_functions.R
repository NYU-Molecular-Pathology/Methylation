runAllBrew <- function(){
    # export RSTUDIO_WHICH_R=/usr/local/bin/R
    system("brew install aspell gdal autoconf automake bash cairo cmake coreutils dos2unix exa fd fontconfig freetype fribidi fzf gawk gcc gdbm gdk-pixbuf gettext giflib git glib gmp gnu-sed gnu-tar gnutls graphite2 harfbuzz htop icu4c ilmbase imagemagick isl jpeg libcroco libde265 libevent libffi libgcrypt libgpg-error libheif libidn2 libmpc libomp libpng librsvg libtasn1 libtermkey libtiff libtool libunistring libuv libvterm libyaml little-cms2 llvm lua lua@5.1 luajit luarocks lzo mpfr msgpack multimarkdown ncurses neovim nettle nspr nss oniguruma openconnect openexr openjpeg openssl@1.1 p11-kit pandoc pandoc-citeproc pandoc-crossref pango pcre pcre2 pdfgrep perl pixman pkg-config poppler qt readline rename ripgrep ruby shared-mime-info shellcheck sqlite stoken swig texinfo tmux tree unbound unibilium vim webp wget x265")
    system('export CFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CCFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CXXFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CPPFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk')
    system("brew install golang libstfl")
    if(dir.exists("/Library/Developer/CommandLineTools")){
        if(length(list.dirs("/Library/Developer/CommandLineTools", recursive = F))>0){
            system("ln -s /Applications/Xcode.app/Contents/Developer /Library/Developer/CommandLineTools")}}
    system("brew unlink gcc && brew link gcc --force")
    system('brew link --overwrite gcc')
    system('brew doctor')
}


fixCompiles <- function(brewExtra=F){
    #system("/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)")
    #cpMakeV()
    system("brew install llvm gcc")
    system("defaults write org.R-project.R NSAppSleepDisabled -bool YES")
    system('echo "export PATH="/usr/local/sbin:$PATH"" >> ~/.zshrc >')
    system('echo "export PATH="/usr/local/sbin:$PATH"" >> ~/.bashrc >')
    system("export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES")
    system('echo "export OBJC_DISABLE_INITIALIZE_FORK_SAFETY="YES"" >> ~/.zshrc')
    system('echo "export OBJC_DISABLE_INITIALIZE_FORK_SAFETY="YES"" >> ~/.bashrc')
    system("xcode-select --install")
    system("xcodebuild -runFirstLaunch")
    system("brew install dpkg openssl libffi go python")
    system('export PATH="/usr/local/opt/python/libexec/bin:$PATH"')
    system('export GOROOT=/usr/local/opt/go/libexec')
    system("export PATH=$PATH:$GOROOT/bin")
    system("brew install cmake gsl gmp libopenmpt open-mpi")
    if (brewExtra == T) {
      runAllBrew()
    }
}
