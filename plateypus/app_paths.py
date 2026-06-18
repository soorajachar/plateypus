#!/usr/bin/env python3

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


APP_NAME = "plateypus"


def running_as_pyinstaller_bundle() -> bool:
    return bool(getattr(sys, "frozen", False)) and hasattr(sys, "_MEIPASS")


def get_resource_directory() -> Path:
    """
    Return the directory containing bundled, read-only plateypus resources.

    In normal pip/editable installs, this is the plateypus package directory.

    In PyInstaller builds, bundled data files are unpacked into sys._MEIPASS.
    If the PyInstaller spec puts package data under a plateypus/ folder, then
    this returns sys._MEIPASS/plateypus.
    """
    if running_as_pyinstaller_bundle():
        bundled_package_dir = Path(sys._MEIPASS) / APP_NAME
        if bundled_package_dir.exists():
            return bundled_package_dir

        return Path(sys._MEIPASS)

    return Path(__file__).resolve().parent


def get_user_data_directory() -> Path:
    """
    Return a user-writable directory for persistent plateypus files.

    Windows:
        C:/Users/<user>/AppData/Roaming/plateypus

    macOS:
        ~/Library/Application Support/plateypus

    Linux:
        ~/.local/share/plateypus
    """
    if os.name == "nt":
        appdata = os.environ.get("APPDATA")
        if appdata:
            user_data_directory = Path(appdata) / APP_NAME
        else:
            user_data_directory = Path.home() / f".{APP_NAME}"
    elif sys.platform == "darwin":
        user_data_directory = Path.home() / "Library" / "Application Support" / APP_NAME
    else:
        xdg_data_home = os.environ.get("XDG_DATA_HOME")
        if xdg_data_home:
            user_data_directory = Path(xdg_data_home) / APP_NAME
        else:
            user_data_directory = Path.home() / ".local" / "share" / APP_NAME

    user_data_directory.mkdir(parents=True, exist_ok=True)
    return user_data_directory


def get_misc_directory() -> Path:
    misc_directory = get_user_data_directory() / "misc"
    misc_directory.mkdir(parents=True, exist_ok=True)
    return misc_directory


def ensure_user_data_files_exist() -> None:
    """
    Create user-writable plateypus misc files if they do not already exist.

    This means users do not need to manually copy, move, or manage anything.
    The app silently initializes itself on first launch.
    """
    resource_directory = get_resource_directory()
    misc_directory = get_misc_directory()

    cytokine_source = resource_directory / "templates" / "cytokineMWDf_template.csv"
    cytokine_destination = misc_directory / "cytokineMWDf.csv"

    kit_source = resource_directory / "templates" / "kitDf_template.csv"
    kit_destination = misc_directory / "kitDf.csv"

    if not cytokine_destination.exists():
        shutil.copyfile(cytokine_source, cytokine_destination)

    if not kit_destination.exists():
        shutil.copyfile(kit_source, kit_destination)
