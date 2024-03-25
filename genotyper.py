#!/usr/bin/env python3
import argparse
import sys
import pysam

# INV genotyper config
__appname__ = "INV genotyper"
__version__ = "0.1"

# INV genotyping classify functions
import classify

def do_work(args):
    # read inv_list.txt to get INV and profile files
    inv_profile = classify.read_inv_list(args.inv_list)

    # open VCF file containing SNP calls
    # get sample info from VCF header
    vcffile=pysam.TabixFile(args.vcf)
    header=vcffile.header
    for line in header:
        if line.split('\t')[0]=="#CHROM":
            sampleinfo=line.split('\t')

    # Read in confident bed file if provided
    if args.bed:
        confident_region = classify.read_bed(args.bed)
    else:
        confident_region = None

    # calculate distance for each INV
    inv_results = {}
    for inv in inv_profile:
        inv_results[inv] = classify.genotype(inv, inv_profile[inv], vcffile, sampleinfo, confident_region, args.skip_missing_tagsnps)

    # write results into output
    classify.write_output(inv_results,args.output, args.min_score)


def parse_arguments(argv):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()

    parser.add_argument("-l", "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    parser.add_argument(
        "--version", action="version", version=__version__
    )

    parser.add_argument(
        "-v", "--vcf", help="Input VCF file containing SNP calls"
    )
    parser.add_argument(
        "-o", "--output", default="inv_genotype.txt",
        help="Output filename (inv_genotype.txt)"
    )

    parser.add_argument(
        "--inv_list", default="inv_list.txt",
        help="File containing a list of INV to genotype"
    )

    parser.add_argument(
        "--min_score", default=0, type=float, help="Minimal confidence score of reporte INV [0]"
    )

    parser.add_argument(
        "--bed", default=None, help="Confident region for VCF file [None]"
    )

    parser.add_argument(
        "--skip_missing_tagsnps", default=False, action='store_true',
        help="Skip missing Tag SNPs [False]"
    )

    return parser.parse_args(argv)


def main(argv):
    """main entry point"""
    # Parse arguments
    args = parse_arguments(argv)

    # Start logging
    #logger.info("%s v%s", pipeline_config.__appname__, pipeline_config.__version__)
    #logger.info(args)
                
    # do work
    do_work(args)

    # done
    #logger.info("Done.")


if __name__ == "__main__":
    main(sys.argv[1:])


