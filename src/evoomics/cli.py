"""Command-line entry point; scientific operations live in the service layer."""
import argparse
from pathlib import Path
import sys

from evoomics import __version__
from evoomics.models import ValidationError
from evoomics.pipeline import run


def main(argv=None):
    parser = argparse.ArgumentParser(prog="evoomics", description="Exploratory evolutionary genomics")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    analyze = commands.add_parser("analyze", help="Aligned DNA FASTA or PASS biallelic SNP VCF to UPGMA")
    analyze.add_argument("--input", type=Path, required=True)
    analyze.add_argument("--format", choices=["fasta", "vcf"], required=True)
    analyze.add_argument("--output", type=Path, required=True)
    analyze.add_argument("--model", choices=["p", "jc69", "k2p"], default="p")
    analyze.add_argument("--heterozygotes", choices=["mask", "iupac"], default="mask")
    analyze.add_argument("--max-site-missing", type=float, default=1.0)
    analyze.add_argument("--min-comparable-sites", type=int, default=1)
    analyze.add_argument("--allow-snp-correction", action="store_true",
                         help="Acknowledge bias; this does NOT correct ascertainment")
    args = parser.parse_args(argv)
    try:
        result = run(args.input, args.output, input_format=args.format, model=args.model,
                     heterozygotes=args.heterozygotes, max_site_missing=args.max_site_missing,
                     min_comparable_sites=args.min_comparable_sites,
                     allow_snp_correction=args.allow_snp_correction)
    except (ValidationError, OSError, UnicodeError, EOFError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"Analysis saved to {result}; read report.md for assumptions and limitations.")
    return 0
