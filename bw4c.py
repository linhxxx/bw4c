#!/usr/bin/env python

import sys
import os
from bw4cModule import convert_reads,as_bam,bam_merge

def main(args=sys.argv[1:]):
    
    if len(args) > 0 and args[0] == "c2t":
        sys.exit(convert_reads(args[1], args[2]))
    elif len(args) > 0 and args[0] == "seqBack":
        # sam input as design for bwa output
        sys.exit(as_bam(args[1]))
    elif len(args) > 0 and args[0] == "bamMerge":
        # sam input as design for bwa output
        sys.exit(bam_merge(args[1], args[2], args[3], args[4]))
    else:
        print("mode c2t fq1 fq2\nPE reads: c2t fq1 fq2\nSE reads: c2t fq1(read1 direction) NA")
        print("mode seqBack sam(or stdin input)")
        print("mode bamMerge watsonBam crickBam outBam, pe_or_se")

if __name__ == "__main__":
    main(sys.argv[1:])
