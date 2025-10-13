import argparse, subprocess, sys, shutil

def main():
    parser = argparse.ArgumentParser(
        description="Run SnakeAltPromoter pipeline (Snakemake wrapper)."
    )
    parser.add_argument("--cores", type=int, default=4, help="CPU cores for Snakemake")
    parser.add_argument("--configfile", type=str, default="config/config.yaml",
                        help="Path to config YAML")
    parser.add_argument("--snakefile", type=str, default="Snakefile",
                        help="Path to Snakefile (default: project root Snakefile)")
    parser.add_argument("--workdir", type=str, default=".",
                        help="Working directory")
    parser.add_argument("--dryrun", action="store_true", help="Snakemake -n")
    args, extra = parser.parse_known_args()

    cmd = [
        shutil.which("snakemake") or "snakemake",
        "--cores", str(args.cores),
        "--snakefile", args.snakefile,
        "--directory", args.workdir,
        "--configfile", args.configfile,
    ]
    if args.dryrun:
        cmd.append("-n")
    cmd.extend(extra)

    print("[sap] Running:", " ".join(cmd))
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)