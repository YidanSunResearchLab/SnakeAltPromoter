import subprocess
import sys


def run_help(module):
    cmd = [sys.executable, "-m", module, "--help"]
    return subprocess.run(cmd, capture_output=True, text=True)


def test_cli_help():
    for module in [
        "snakealtpromoter.workflows.Snakealtpromoter",
        "snakealtpromoter.workflows.Genomesetup",
        "snakealtpromoter.cli",
        "snakealtpromoter.ui.launch",
    ]:
        result = run_help(module)
        assert result.returncode == 0, result.stderr
