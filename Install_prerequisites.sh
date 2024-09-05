#!/bin/bash

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to install prerequisites on Linux (Debian-based)
install_linux_debian() {
    echo "Updating package list..."
    sudo apt-get update

    echo "Installing prerequisites..."
    sudo apt-get install -y curl make build-essential

    echo "Prerequisites installed on Debian-based Linux."
}

# Function to install prerequisites on Linux (Red Hat-based)
install_linux_redhat() {
    echo "Updating package list..."
    sudo yum check-update

    echo "Installing prerequisites..."
    sudo yum install -y curl make gcc gcc-c++

    echo "Prerequisites installed on Red Hat-based Linux."
}

# Function to install prerequisites on macOS
install_macos() {
    if ! command_exists brew; then
        echo "Homebrew not found. Installing Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    fi

    echo "Installing prerequisites..."
    brew install curl make gcc

    echo "Prerequisites installed on macOS."
}

# Function to guide Windows users
install_windows() {
    echo "For Windows, please install the following manually:"
    echo "1. Download and install Git for Windows: https://gitforwindows.org/"
    echo "2. Download and install Chocolatey (a package manager for Windows): https://chocolatey.org/install"
    echo "3. Use Chocolatey to install prerequisites by running: choco install curl make gcc"
}

# Determine the OS and install prerequisites accordingly
case "$(uname -s)" in
    Linux)
        # Determine if Debian-based or Red Hat-based
        if [ -f /etc/debian_version ]; then
            install_linux_debian
        elif [ -f /etc/redhat-release ]; then
            install_linux_redhat
        else
            echo "Unsupported Linux distribution."
            exit 1
        fi
        ;;
    Darwin)
        install_macos
        ;;
    CYGWIN*|MINGW32*|MINGW64*|MSYS*|MINGW*)
        install_windows
        ;;
    *)
        echo "Unsupported OS."
        exit 1
        ;;
esac
