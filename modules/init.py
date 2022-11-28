import os
import argparse

program_name = "VJ-gene-finder"
program_version = "1.0"
program_description = ("Extraction of chicken TCR VJ gene segments"
                       + "from a reference genome")
program_epilog = "https://github.com/simonfrueh/VJ-gene-finder"


# Parse command line arguments
def parse_arguments():
    # Initialize command line argument parser
    parser = argparse.ArgumentParser(
        prog=program_name,
        description=program_description,
        epilog=program_epilog)

    # Add and parse arguments
    parser.add_argument("filename", help="Input file")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s "+program_version)
    parser.add_argument("-o", "--output", help="Output directory")
    args = parser.parse_args()

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

    return path, filename, output_dir
