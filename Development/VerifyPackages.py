#!/usr/bin/env python3
# Script name: VerifyPackages.py
# Purpose: Ensures pip and python modules are loaded and installed
# Date Created: February 28, 2025
# Version: 1.0
# Author: Jonathan Serrano
# Copyright (c) NYULH Jonathan Serrano, 2025

import subprocess
import sys
import logging
import argparse
from pathlib import Path

def create_virtual_environment(venv_path: Path) -> None:
    """Create a virtual environment if it doesn't exist."""
    if not venv_path.exists():
        logging.info(f"Creating virtual environment at {venv_path}")
        subprocess.run([sys.executable, "-m", "venv", str(venv_path)], check=True)
    else:
        logging.info(f"Virtual environment already exists at {venv_path}")

def upgrade_pip(venv_path: Path) -> None:
    """Upgrade pip to the latest version."""
    logging.info("Upgrading pip to the latest version.")
    pip_executable = venv_path / ("Scripts" if sys.platform == "win32" else "bin") / "pip"
    subprocess.run([str(pip_executable), "install", "--upgrade", "pip"], check=True)

def install_package(venv_path: Path, package: str) -> None:
    """Install a package using pip within the virtual env."""
    logging.info(f"Installing '{package}' in the virtual environment.")
    pip_executable = venv_path / ("Scripts" if sys.platform == "win32" else "bin") / "pip"
    subprocess.run([str(pip_executable), "install", package], check=True)

def parse_arguments():
    """Parse args to get list of packages to install & venv path."""
    parser = argparse.ArgumentParser(
        description="Install required Python packages in a virtual environment."
    )
    parser.add_argument("packages", nargs="+", help="List of packages to install.")
    parser.add_argument(
        "--venv", type=Path, default=Path("./venv"),
        help="Path to the virtual environment directory (default: ./venv)"
    )
    return parser.parse_args()

def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    args = parse_arguments()
    venv_path = args.venv

    create_virtual_environment(venv_path)
    upgrade_pip(venv_path)

    for package in args.packages:
        try:
            install_package(venv_path, package)
        except subprocess.CalledProcessError as error:
            logging.critical(f"Failed to install '{package}': {error}")
            sys.exit(1)
    print("Be sure to include the following to your python script:")
    print("#!/usr/bin/env ./venv/bin/python")

if __name__ == "__main__":
    main()
