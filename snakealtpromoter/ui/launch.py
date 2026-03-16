#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import subprocess
from importlib import resources


def main():
    if any(arg in ("-h", "--help") for arg in sys.argv[1:]):
        print("usage: sap-ui")
        print("Launch the SnakeAltPromoter Streamlit UI.")
        return

    # locate the app file inside the installed package
    with resources.path("snakealtpromoter.ui", "stl_app.py") as app_path:
        cmd = [
            "streamlit",
            "run",
            str(app_path),
            "--browser.gatherUsageStats",
            "false"
        ]
        subprocess.run(cmd)


if __name__ == "__main__":
    main()
