import subprocess
import sys
import pkg_resources


def install_pip(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])


def install_conda(package):
    subprocess.check_call(["conda", "install", "-y", package])


def check_and_install(package):
    try:
        pkg_resources.get_distribution(package)
    except pkg_resources.DistributionNotFound:
        try:
            install_pip(package)
        except subprocess.CalledProcessError:
            print(f"Error: Failed to install '{package}' using pip. Trying with conda...")
            try:
                install_conda(package)
            except subprocess.CalledProcessError:
                print(f"Error: Failed to install '{package}' using conda. Please install it manually.")
                sys.exit(1)


# Example usage:
packages_to_install = ["methylprep", "methylcheck"]

for package in packages_to_install:
    check_and_install(package)

