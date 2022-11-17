import os
import argparse


# Parse command line arguments
def parse_arguments():
    # Initialize command line argument parser
    parser = argparse.ArgumentParser(
        prog="VJ-gene-finder",
        description="Extraction of chicken TCR VJ gene segments"
                    + "from a reference genome",
        epilog="https://github.com/simonfrueh/VJ-gene-finder")

    # Add and parse arguments
    # To Do:
    # - option for output directory
    # - option for version information
    parser.add_argument("filename", help="Input file")
    args = parser.parse_args()

    if args.filename:
        print("Input file: %s" % args.filename)

    path = os.path.dirname(args.filename).replace(" ", "")

    # Add trailing slash for later usage
    if len(path) > 0:
        if path[-1] != "/":
            path += "/"

    filename = os.path.basename(args.filename)

    return path, filename
