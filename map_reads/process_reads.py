__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "processes unmapped reads"

import subprocess

def _check_subprocess_run(returncode, stderrdata, runinfo):
    if returncode == 0:
        print("STATUS: "+runinfo+" complete")
    else:
        print("WARNING: Issue with "+runinfo+" reads")
        print(stderrdata)

def process_unmapped_reads(args, samtoolspath):
    try: 
        os.mkdir(args.output_path+"_tmp")
    except:
        pass
    # -- pull out unmapped reads
    args=[samtoolspath+"samtools", "view", "-@", args.processors, "-b", "-h", "-f", "4", args.path2possorted_genome_bam, "-o",  args.output_path+"_tmp"+"/unmmapped.bam"]
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    print("STATUS: Extracting unmapped reads..")
    stdoutdata, stderrdata = process.communicate()
    _check_subprocess_run(process.returncode, stderrdata, "extracting unmapped")
    
    # -- convert to fastq file
    args = [samtoolspath+"samtools", "bam2fq", "-n","-O", "-s", args.output_path+"_tmp"+"ummapped.fq", args.output_path+"_tmp"+"unmmapped.bam"]
    f = open(args.output_path+"_tmp"+"ummapped.fq" "w") 
    process = subprocess.Popen(args, stdout=f, stderr=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    f.close()
    _check_subprocess_run(process.returncode, stderrdata, "converting reads to fastq")