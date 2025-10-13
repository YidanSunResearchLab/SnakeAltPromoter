# streamlit_app.py
import streamlit as st
import subprocess, sys, os, shlex

st.set_page_config(page_title="SnakeAltPromoter UI", layout="wide")
st.title("SnakeAltPromoter — Local UI")

st.caption(
    "UI interface for genome setup and main pipeline. "
)

tab1, tab2 = st.tabs(["Genome setup (Genomesetup)", "Main pipeline (Snakealtpromoter.py)"])

def run_and_stream(cmd_list):
    st.write("**Command to run:**")
    st.code(" ".join(shlex.quote(c) for c in cmd_list))
    with subprocess.Popen(
        cmd_list,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1
    ) as p:
        log = st.empty()
        lines = []
        for line in p.stdout:
            lines.append(line.rstrip())
            # show the last 300 lines to keep the UI responsive
            log.code("\n".join(lines[-300:]))
        ret = p.wait()
    if ret == 0:
        st.success("Done.")
    else:
        st.error(f"Failed (exit code {ret}).")

with tab1:
    st.subheader("Genome setup — build genome index and annotations")
    organism = st.text_input("organism", value="hg38")
    fasta = st.text_input("organism_fasta (path to genome.fa)", value="/path/to/your/genome.fa")
    gtf = st.text_input("genes_gtf (path to genes.gtf)", value="/path/to/your/genes.gtf")
    out_dir = st.text_input("Output directory (-o)", value="./genome")
    threads = st.number_input("Threads (--threads)", min_value=1, max_value=128, value=30)
    if st.button("Start build", type="primary"):
        cmd = [
            "Genomesetup",
            "--organism", organism,
            "--organism_fasta", fasta,
            "--genes_gtf", gtf,
            "-o", out_dir,
            "--threads", str(threads),
        ]
        run_and_stream(cmd)

with tab2:
    st.subheader("Snakealtpromoter.py — main RNA-seq pipeline")
    input_dir = st.text_input("-i Input FASTQ directory", value="/path/to/input/fastqs/dir/")
    genome_dir = st.text_input("--genome_dir (output from the previous step)", value="/path/to/genomesetup/dir/")
    out_dir2 = st.text_input("-o Output directory", value="/path/to/output/dir/")
    threads2 = st.number_input("--threads", min_value=1, max_value=128, value=30)
    organism2 = st.text_input("--organism", value="hg38")
    ctrl = st.text_input("--control_condition", value="Healthy")
    case = st.text_input("--test_condition", value="Heart_Failure")

    extra = st.text_input(
        "Optional: extra advanced args to pass through verbatim",
        value=""
    )
    if st.button("Run pipeline", type="primary"):
        cmd = [
            "Snakealtpromoter.py",
            "-i", input_dir,
            "--genome_dir", genome_dir,
            "-o", out_dir2,
            "--threads", str(threads2),
            "--organism", organism2,
            "--control_condition", ctrl,
            "--test_condition", case,
        ]
        if extra.strip():
            cmd += shlex.split(extra)
        run_and_stream(cmd)