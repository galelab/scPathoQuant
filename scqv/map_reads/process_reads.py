__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "processes unmapped reads"

import re
import os
import subprocess
from sys import platform
import pysam
PATH = os.path.dirname(os.path.abspath(__file__))
PATH = re.sub("map_reads", "", PATH)
if platform == "linux":
    bowtie2path = os.path.join(PATH, "aligntools", "bowtie2-2.4.2-linux-x86_64")
elif platform == "OS":
    bowtie2path = os.path.join(PATH, "aligntools","bowtie2-2.4.2-macos-x86_64")
else:
    ValueError("Program wont run on this operating system "+platform)

def _check_subprocess_run(returncode, stderrdata, runinfo):
    if returncode == 0:
        print("STATUS: "+runinfo+" complete")
    else:
        print("WARNING: Issue with "+runinfo+" reads")
        print(stderrdata)

def process_unmapped_reads(args, viable_cb):
    try: 
       os.mkdir(args.output_path)
    except:
        pass 
    try: 
        print("STATUS: Generating _tmp/ files")
        os.mkdir(os.path.join(args.output_path, "_tmp"))
        # -- pull out unmapped reads

        pysam.view( "-@", str(args.processors), "-b", "-h", "-f", "4", os.path.join(args.path10x, "outs", "possorted_genome_bam.bam"), "-o", os.path.join(args.output_path, "_tmp","unmapped.bam"), catch_stdout=False)
        pysam.view("-@", str(args.processors), "-h", os.path.join(args.output_path,"_tmp", "unmapped.bam"), "-o", os.path.join(args.output_path, "_tmp", "unmapped.sam"), catch_stdout=False)

        # arg=[ os.path.join(samtoolspath, "samtools"), "view", "-@", args.processors,
        #     "-b", "-h", "-f", "4", os.path.join(args.path10x, "outs", "possorted_genome_bam.bam"),
        #     "-o",  os.path.join(args.output_path, "_tmp","unmapped.bam")]
        # process = subprocess.Popen(arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        # print("STATUS: Extracting unmapped reads..")
        # stdoutdata, stderrdata = process.communicate()
        # _check_subprocess_run(process.returncode, stderrdata, "extracting unmapped")
        
        # -- convert unmapped bam to sam file
        # arg=[ os.path.join(samtoolspath, "samtools"),"view", "-@", args.processors, 
        #     "-h", os.path.join(args.output_path,"_tmp", "unmapped.bam"),
        #     "-o", os.path.join(args.output_path, "_tmp", "unmapped.sam")]
        # process = subprocess.Popen(arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        # print("STATUS: Convert bam to sam for unmmaped reads..")
        # stdoutdata, stderrdata = process.communicate()
        # _check_subprocess_run(process.returncode, stderrdata, "convert bam to sam")

        # -- pull out umi and cell barcode information 
        print ("STATUS: extracting cell barcodes and UMIs...")
        with open(os.path.join(args.output_path, "_tmp", "barcode_umi_read_table.csv"), "w") as fout:
            fout.write("cell_barcode,umi,read\n")
            with open(os.path.join(args.output_path,"_tmp","unmapped.sam")) as fin:
                for line in fin:
                    if line.startswith("@"):
                        pass
                    else: 
                        line = line.strip()
                        larray = line.split("\t")
                        cellbarcode=False
                        umi=False
                        for i in larray:
                            if i.startswith("CB"):
                                cellbarcode=i
                            elif i.startswith("UB"):
                                umi=i
                        if cellbarcode != False and umi != False and re.sub("CB:Z:", "", cellbarcode) in viable_cb: 
                            fout.write(re.sub("CB:Z:", "", cellbarcode)+","+re.sub("UB:Z:", "", umi)+","+larray[0]+"\n")

        # -- convert to fastq file
        pysam.bam2fq('-n', '-0', os.path.join(args.output_path, "_tmp", "unmapped.fq"), os.path.join(args.output_path,"_tmp","unmapped.bam"),  catch_stdout=False)
        # arg = [os.path.join(samtoolspath, "samtools"), "bam2fq", "-@", args.processors, "-n","-O", "-s", os.path.join(args.output_path, "_tmp", "unmapped.fq"), os.path.join(args.output_path,"_tmp","unmapped.bam")]
        # f = open(os.path.join(args.output_path, "_tmp", "unmapped.fq"), "w") 
        # process = subprocess.Popen(arg, stdout=f, stderr=subprocess.PIPE)
        # stdoutdata, stderrdata = process.communicate()
        # f.close()
        # _check_subprocess_run(process.returncode, stderrdata, "converting reads to fastq")
    except:
        print("STATUS: _tmp folder already exists so not regenerating data, path "+str(os.path.join(args.output_path, "_tmp")))
        pass