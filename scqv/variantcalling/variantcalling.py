__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "map reads with bowtie (currently only version)"

import re
import os
import glob
import subprocess
import scqv.extra_functions as ef
from sys import platform


PATH = os.path.dirname(os.path.abspath(__file__))
PATH = re.sub("variantcalling", "", PATH)
if platform == "linux":
    samtoolspath = os.path.join(PATH, "extra_tools", "samtoolsv1.11_linux","bin")
    bcftoolspath =  os.path.join(PATH, "extra_tools", "bcftools1.10.2_linux", "bin")
    vcfutilspath =  os.path.join(PATH, "extra_tools", "vcfutils", "bin")


def variantcaller(args):
    files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
    genomefile = files_genome[0]

    arg=[os.path.join(samtoolspath, "samtools"), "faidx", genomefile]
    ef._run_subprocesses(arg, "STATUS: indexing genome file", "indexing genome file")

    arg=[os.path.join(samtoolspath, "samtools"), "mpileup", "-f ", genomefile, os.path.join(args.output_path,"virus_al_sort.bam"), "-o", os.path.join(args.output_path,"virus.bcf")]
    ef._run_subprocesses(arg, "STATUS: getting snps", "getting snps")

    arg=[os.path.join(bcftoolspath, "bcftools"), "view",  os.path.join(args.output_path,"virus.bcf"), "-o", os.path.join(args.output_path,"virus_var.bcf")]
    ef._run_subprocesses(arg, "STATUS: Calling snps", "calling snps")

    arg=[os.path.join(bcftoolspath, "bcftools"), "view",  os.path.join(args.output_path,"virus_var.bcf"), "|", os.path.join(vcfutilspath, "vcfutils.pl"), "varFilter", "-d", 2, ">",  os.path.join(args.output_path,"virus_var_final.vcf")]
    ef._run_subprocesses(arg, "STATUS: Converting to vcf", "converting to vcf")
