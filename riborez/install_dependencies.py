#!/usr/bin/env python3
"""
Dependency installation for RiboRez
"""

import os
import sys
import subprocess
import platform
import urllib.request
import shutil
from pathlib import Path

def check_datasets_installed():
    """Check if datasets command is available."""
    try:
        result = subprocess.run(['datasets', '--version'], 
                              capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def install_via_conda():
    """Try installing via conda."""
    # Check if conda is available
    if not shutil.which('conda'):
        print("[INFO] Conda not found, skipping conda installation")
        return False
    
    print("[INFO] Attempting conda installation...")
    try:
        result = subprocess.run(
            ["conda", "install", "-c", "conda-forge", "ncbi-datasets-cli", "-y"],
            capture_output=True, text=True, check=True
        )
        print("[SUCCESS] Conda installation completed")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[FAILED] Conda installation failed: {e.stderr}")
        return False
    except Exception as e:
        print(f"[ERROR] Conda installation error: {e}")
        return False

def install_via_download():
    """Install NCBI Datasets CLI binary via direct download."""
    print("[INFO] Attempting direct download installation...")
    
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
    try:
        datasets_path = bin_dir / datasets_name
        urllib.request.urlretrieve(datasets_url, datasets_path)
        datasets_path.chmod(0o755)  # Make executable
        print(f"[SUCCESS] Downloaded datasets to {datasets_path}")
    except Exception as e:
        print(f"[ERROR] Failed to download datasets: {e}")
        return False
    
    # Download dataformat
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
        print(f"[INFO] Please add the following to your shell profile (.bashrc, .zshrc, etc.):")
        print(f"export PATH=\"$HOME/.local/bin:$PATH\"")
        # Also add to current session
        os.environ['PATH'] = f"{bin_dir}:{path_env}"
    
    print("[SUCCESS] Direct download installation completed!")
    return True

def install_ncbi_datasets():
    """Install NCBI Datasets CLI using multiple methods."""
    if check_datasets_installed():
        return True
    
    print("[INFO] NCBI Datasets CLI not found. Attempting installation...")
    
    # Try conda first (more reliable)
    if install_via_conda():
        if check_datasets_installed():
            return True
        else:
            print("[WARNING] Conda installation completed but datasets command not found in PATH")
    
    # Fall back to direct download
    if install_via_download():
        if check_datasets_installed():
            return True
        else:
            print("[WARNING] Direct download completed but datasets command not found in PATH")
    
    return False

def ensure_dependencies():
    """Ensure all dependencies are installed."""
    if not check_datasets_installed():
        return install_ncbi_datasets()
    return True 