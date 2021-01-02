__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "processes unmapped reads"

import re
import os
import subprocess

def _check_subprocess_run(returncode, stderrdata, runinfo):
    if returncode == 0:
        print("STATUS: "+runinfo+" complete")
    else:
        print("WARNING: Issue with "+runinfo+" reads")
        print(stderrdata)

def process_unmapped_reads(args, samtoolspath, viable_cb):
    try: 
        os.mkdir(args.output_path+"_tmp")
    except:
        pass
    # -- pull out unmapped reads
    args=[samtoolspath+"samtools", "view", "-@", args.processors, "-b", "-h", "-f", "4", os.path.join(args.path10x, "outs", "possorted_genome_bam.bam"), "-o",  os.path.join(args.output_path, "_tmp","unmmapped.bam")]
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    print("STATUS: Extracting unmapped reads..")
    stdoutdata, stderrdata = process.communicate()
    _check_subprocess_run(process.returncode, stderrdata, "extracting unmapped")
    
    # -- convert unmapped bam to sam file
    args=[samtoolspath+"samtools","view" "-@", args.processors, "-h", os.path.join(args.output_path,"_tmp", "unmapped.bam"), "-o", os.path.join(args.output_path, "_tmp", "unmapped.sam")]
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    print("STATUS: Convert bam to sam for unmmaped reads..")
    stdoutdata, stderrdata = process.communicate()
    _check_subprocess_run(process.returncode, stderrdata, "convert bam to sam")

    # -- pull out umi and cell barcode information 
    print ("STATUS: extracting cell barcodes and UMIs...")
    with open(os.path.join(args.output_path, "_tmp"," barcode_umi_read_table.csv"), "w") as fout:
        fout.write("cell_barcode,umi,read\n")
        with open(args.output_path+"_tmp"+"/unmapped.sam") as fin:
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
    args = [samtoolspath+"samtools", "bam2fq", "-@", args.processors, "-n","-O", "-s", os.path.join(args.output_path, "_tmp", "ummapped.fq"), os.path.join(args.output_path,"_tmp","unmmapped.bam")]
    f = open(args.output_path+"_tmp"+"ummapped.fq" "w") 
    process = subprocess.Popen(args, stdout=f, stderr=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    f.close()
    _check_subprocess_run(process.returncode, stderrdata, "converting reads to fastq")