#!/usr/bin/env python3
"""
Installation script for NCBI Datasets CLI
Tries multiple installation methods to ensure the tool is available.
"""

import os
import sys
import subprocess
import platform
import urllib.request
import shutil
from pathlib import Path

def run_command(cmd, description):
    """Run a command and return success status."""
    print(f"[INFO] {description}...")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"[SUCCESS] {description}")
            return True
        else:
            print(f"[FAILED] {description}: {result.stderr}")
            return False
    except Exception as e:
        print(f"[ERROR] {description}: {e}")
        return False

def check_datasets_installed():
    """Check if datasets command is available."""
    try:
        result = subprocess.run(['datasets', '--version'], 
                              capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def install_via_pip():
    """Try installing via pip."""
    return run_command("pip install ncbi-datasets-cli", "Installing via pip")

def install_via_conda():
    """Try installing via conda."""
    # Check if conda is available
    if not shutil.which('conda'):
        print("[INFO] Conda not found, skipping conda installation")
        return False
    
    # Try to install via conda-forge
    return run_command("conda install -c conda-forge ncbi-datasets-cli -y", 
                      "Installing via conda")

def download_binary():
    """Download the binary directly from NCBI."""
    system = platform.system().lower()
    machine = platform.machine().lower()
    
    # Determine the correct URL
    if system == "linux" and "x86_64" in machine:
        datasets_url = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets"
        dataformat_url = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat"
        datasets_name = "datasets"
        dataformat_name = "dataformat"
    elif system == "darwin":  # macOS
        datasets_url = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets"
        dataformat_url = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat"
        datasets_name = "datasets"
        dataformat_name = "dataformat"
    elif system == "windows":
        datasets_url = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/datasets.exe"
        dataformat_url = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/dataformat.exe"
        datasets_name = "datasets.exe"
        dataformat_name = "dataformat.exe"
    else:
        print(f"[WARNING] Unsupported platform: {system} {machine}")
        return False
    
    # Create bin directory in user's home
    bin_dir = Path.home() / ".local" / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    
    # Download datasets
    print(f"[INFO] Downloading datasets binary...")
    try:
        datasets_path = bin_dir / datasets_name
        urllib.request.urlretrieve(datasets_url, datasets_path)
        datasets_path.chmod(0o755)  # Make executable
        print(f"[SUCCESS] Downloaded datasets to {datasets_path}")
    except Exception as e:
        print(f"[ERROR] Failed to download datasets: {e}")
        return False
    
    # Download dataformat
    print(f"[INFO] Downloading dataformat binary...")
    try:
        dataformat_path = bin_dir / dataformat_name
        urllib.request.urlretrieve(dataformat_url, dataformat_path)
        dataformat_path.chmod(0o755)  # Make executable
        print(f"[SUCCESS] Downloaded dataformat to {dataformat_path}")
    except Exception as e:
        print(f"[ERROR] Failed to download dataformat: {e}")
        return False
    
    # Add to PATH if not already there
    path_env = os.environ.get('PATH', '')
    if str(bin_dir) not in path_env:
        print(f"[INFO] Adding {bin_dir} to PATH")
        print(f"[INFO] Please add the following to your shell profile (.bashrc, .zshrc, etc.):")
        print(f"export PATH=\"$HOME/.local/bin:$PATH\"")
    
    return True

def main():
    """Main installation function."""
    print("=== NCBI Datasets CLI Installation ===")
    
    # Check if already installed
    if check_datasets_installed():
        print("[SUCCESS] NCBI Datasets CLI is already installed!")
        return True
    
    print("[INFO] NCBI Datasets CLI not found. Attempting installation...")
    
    # Try different installation methods
    methods = [
        ("pip", install_via_pip),
        ("conda", install_via_conda),
        ("binary download", download_binary)
    ]
    
    for method_name, method_func in methods:
        print(f"\n--- Trying {method_name} ---")
        if method_func():
            # Verify installation
            if check_datasets_installed():
                print(f"\n[SUCCESS] NCBI Datasets CLI installed successfully via {method_name}!")
                return True
            else:
                print(f"[WARNING] Installation via {method_name} completed but datasets command not found")
    
    print("\n[ERROR] All installation methods failed!")
    print("\nManual installation options:")
    print("1. Install via conda: conda install -c conda-forge ncbi-datasets-cli")
    print("2. Download manually from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/")
    print("3. Use the download script provided in the documentation")
    
    return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 