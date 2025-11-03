# import subprocess, sys

# def main():
#     app_module = "app/stl_app.py"
#     cmd = ["streamlit", "run", app_module, "--browser.gatherUsageStats", "false"]
#     sys.exit(subprocess.call(cmd))

# ui/launch.py
from importlib.resources import files, as_file
import sys

def main():
    # import here so missing deps give a clean error
    import streamlit.web.cli as stcli

    # locate the app file inside the installed package
    app_res = files("ui").joinpath("stl_app.py")

    # if package is zipped, as_file() creates a real temp path for Streamlit
    with as_file(app_res) as app_path:
        sys.argv = ["streamlit", "run", str(app_path), "--browser.gatherUsageStats", "false"]
        raise SystemExit(stcli.main())
