__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "map reads with bowtie (currently only version)"

import glob
import subprocess



def _check_subprocess_run(returncode, stderrdata, runinfo):
    if returncode == 0:
        print("STATUS: "+runinfo+" complete")
    else:
        print("WARNING: Issue with "+runinfo+" reads")
        print(stderrdata)

def _run_subprocesses(args, status, error_message):
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    print("STATUS: generating libraries ..")
    stdoutdata, stderrdata = process.communicate()
    _check_subprocess_run(process.returncode, stderrdata, "extracting unmapped")

def map2viralgenome(args, bowtie2path, samtoolspath):

    # -- check if libraries are present if not generate 
    files = glob.glob(args.path2genome+"/*.bt2")
    if len(files) == 0:
        files_genome = glob.glob(args.path2genome+".fa")
        if len(files_genome) > 1:
            print ("WARNING:  two fasta files in genome folder.")
        arg=[bowtie2path+"bowtie2-build", files_genome[0], args.path2genome+"genome"]
        _run_subprocesses(arg, "STATUS: generating libraries ...", "extracting unmapped")

    else:
        print("STATUS: libraries have already been made for this genome")

    # -- align reads 
    arg=[bowtie2path+"bowtie2", "-q", args.output_path+"_tmp/unmmapped.fq", "-x", args.path2genome+"genome", "-S", args.output_path+"virus_aligned.sam"]
    _run_subprocesses(arg, "STATUS: Align reads ...", "aligning reads")
    
    # --generate necessary output files
    arg=[samtoolspath+"samtools", "view", "-F", "4", "-Sb", args.output_path+"virus_aligned.sam", "-o", args.output_path+"virus_aligned.bam"]
    _run_subprocesses(arg, "STATUS: Extracting mapped reads and converting to bam", "extracting mapped reads and converting to bam")

    arg=[samtoolspath+"samtools", "sort", args.output_path+"virus_aligned.bam", "-o", args.output_path+"viral_aligned_sort.bam"]
    _run_subprocesses(arg, "STATUS: Sorting bam", "sorting bam")

    arg=[samtoolspath+"samtools", "index", args.output_path+"virus_aligned_sort.bam", "viral_aligned_sort.bam.bai"]
    _run_subprocesses(arg, "STATUS: generating bam index", "generating bam index")


