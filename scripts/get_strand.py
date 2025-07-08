#!/usr/bin/env python3

import sys

def get_strand(strand_file):
    with open(strand_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "Fraction of reads explained by" in line:
                if "++,--" in line and float(line.split(":")[1].strip()) > 0.7:
                    return "1"  # Forward stranded
                elif "+-,-+" in line and float(line.split(":")[1].strip()) > 0.7:
                    return "2"  # Reverse stranded
    return "0"  # Unstranded

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python get_strand.py <strand.txt>")
    print(get_strand(sys.argv[1]))
