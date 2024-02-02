import os
import argparse

program_name = "VJ-gene-finder"
program_version = "0.1+"
program_description = ("Extraction of chicken TCR VJ gene segments"
                       + "from a reference genome")
program_epilog = "https://github.com/simonfrueh/VJ-gene-finder"


# Parse command line arguments
def parse_arguments():
    # Initialize command line argument parser
    p = argparse.ArgumentParser(prog=program_name,
                                description=program_description,
                                epilog=program_epilog)

    # Add and parse arguments
    p.add_argument("filename", help="Input file")
    p.add_argument("-v", "--version", action="version",
                   version="%(prog)s "+program_version)
    p.add_argument("-o", "--output", help="Output directory")
    p.add_argument("-st",
                   "--skip_trgf",
                   action="store_true",
                   help="Skip assignment of TR family name (based on amino "
                   "acid motif) for V genes and use default names instead "
                   "(TRV-SEL for V genes with single-exon leader peptide,"
                   " TRV for V genes with two-exon leader peptide).")
    p.add_argument("-ss",
                   "--skip_selp",
                   action="store_true",
                   help="Skip search for V genes with single-exon leader "
                   "peptide.")
    args = p.parse_args()

    # Process Input filename (mandatory)
    print("Input file: %s" % args.filename)
    path = os.path.dirname(args.filename)

    # Add trailing slash for later usage
    if len(path) > 0:
        if path[-1] != "/":
            path += "/"
    filename = os.path.basename(args.filename)

    output_dir = path + "output/"
    if args.output:
        output_dir = args.output
        print("Output directory: %s" % output_dir)

    # Add trailing slash for later usage
    if output_dir[-1] != "/":
        output_dir += "/"

    return path, filename, output_dir, args.skip_trgf, args.skip_selp
