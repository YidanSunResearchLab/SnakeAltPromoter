#!/usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
from importlib import resources


def main():
    # locate the app file inside the installed package
    with resources.path('snakealtpromoter.ui', 'stl_app.py') as app_path:
        cmd = [
            "streamlit",
            "run",
            str(app_path),
            "--browser.gatherUsageStats",
            "false"
        ]
        subprocess.run(cmd)


if __name__ == '__main__':
    main()
