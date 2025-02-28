#!/usr/bin/env python3
# Script name: VerifyPackages.py
# Purpose: Ensures pip and python modules are loaded and installed
# Date Created: February 28, 2025
# Version: 1.0
# Author: Jonathan Serrano
# Copyright (c) NYULH Jonathan Serrano, 2025


import subprocess
import sys
import pkg_resources
import logging
import shutil
import platform
import argparse
import importlib.util
from typing import List


def ensure_and_update_pip() -> None:
    """
    Ensure that pip is installed; if not, bootstrap pip using ensurepip.
    Then, upgrade pip to the latest version.
    """
    try:
        # Check if pip is available
        import pip  # noqa: F401
        logging.info("pip is available.")
    except ImportError:
        logging.info("pip not found. Bootstrapping pip using ensurepip.")
        try:
            import ensurepip
            ensurepip.bootstrap()
            logging.info("pip bootstrapped successfully.")
        except Exception as bootstrap_error:
            logging.critical(f"Failed to bootstrap pip: {bootstrap_error}")
            sys.exit(1)

    # Upgrade pip to the latest version using the current Python interpreter
    pip_upgrade_cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "pip"]
    logging.info("Upgrading pip to the latest version...")
    try:
        subprocess.run(pip_upgrade_cmd, check=True)
        logging.info("pip upgraded successfully.")
    except subprocess.CalledProcessError as upgrade_error:
        logging.error(f"Failed to upgrade pip: {upgrade_error}")
        # Continue execution even if pip upgrade fails


def get_pip_command() -> List[str]:
    """
    Return the pip command list that ensures the pip executable corresponds
    to the current Python interpreter.
    """
    system = platform.system()
    logging.info(f"Operating System detected: {system}")
    return [sys.executable, "-m", "pip"]


def install_package_with_pip(package: str) -> None:
    """
    Install a package using pip.
    """
    pip_cmd = get_pip_command() + ["install", package]
    logging.info(f"Installing '{package}' via pip with command: {' '.join(pip_cmd)}")
    subprocess.run(pip_cmd, check=True)


def install_package_with_conda(package: str) -> None:
    """
    Install a package using conda.
    """
    if shutil.which("conda") is None:
        raise EnvironmentError("Conda is not available on this system.")
    logging.info(f"Installing '{package}' via conda.")
    subprocess.run(["conda", "install", "-y", package], check=True)


def check_and_install(package: str) -> None:
    """
    Check if a package is installed. If not, attempt installation via pip,
    and then try conda if pip fails.
    """
    try:
        pkg_resources.get_distribution(package)
        logging.info(f"Package '{package}' is already installed.")
    except pkg_resources.DistributionNotFound:
        logging.info(f"Package '{package}' not found. Attempting installation via pip.")
        try:
            install_package_with_pip(package)
        except subprocess.CalledProcessError as pip_error:
            logging.error(f"Pip installation failed for '{package}': {pip_error}. Trying conda...")
            try:
                install_package_with_conda(package)
            except (subprocess.CalledProcessError, EnvironmentError) as conda_error:
                logging.error(f"Conda installation also failed for '{package}': {conda_error}.")
                raise RuntimeError(f"Failed to install '{package}'. Please install it manually.") from conda_error


def parse_arguments() -> List[str]:
    """
    Parse command-line arguments to obtain a list of packages to install.
    """
    parser = argparse.ArgumentParser(
        description="Install required Python packages using pip and conda."
    )
    parser.add_argument("packages", nargs="+", help="List of packages to install.")
    args = parser.parse_args()
    return args.packages


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    
    # Ensure pip is installed and updated before proceeding
    ensure_and_update_pip()

    packages_to_install = parse_arguments()

    for package in packages_to_install:
        try:
            check_and_install(package)
        except Exception as error:
            logging.critical(error)
            sys.exit(1)


if __name__ == "__main__":
    main()
