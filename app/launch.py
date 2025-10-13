import subprocess, sys

def main():
    app_module = "app/stl_app.py"
    cmd = ["streamlit", "run", app_module, "--browser.gatherUsageStats", "false"]
    sys.exit(subprocess.call(cmd))