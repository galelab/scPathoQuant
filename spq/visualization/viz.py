__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "visualize counts"

import os
import glob
import re
import subprocess
import pysam
import pandas as pd
import seaborn as sns
import spq.extra_functions as ef
from spq.map_reads import map_reads as mr
from spq.quantify import patho_copies as vc

def process_ini_file(args, filename):
    with open(filename, "r") as fin:
        outputfile = re.sub(".ini$", "edit.ini", filename)
        with open(outputfile, 'w') as fout:
            for line in fin:
                if line.startswith("labels = false"):
                    fout.write("labels = true\n")
                elif line.startswith("color = #666666"):
                    fout.write("color = green\n")
                else:
                    fout.write(line)

def generate_viral_copy_plots(args, dfumidict):
    print ("STATUS: Generating violin plots for viral copies")
    for virus_name, dfumi in dfumidict.items():
        if dfumi.shape[0] > 0:
            virusname_rename = re.sub("\/", "_", virus_name )
            virusname_rename = re.sub(" ", "_",virusname_rename)
            sns.set_theme(style="whitegrid")
            dfumi[virus_name]=[virus_name]*dfumi.shape[0]
            ax = sns.violinplot(x=virus_name, y="umi",data=dfumi)
            ax.set(xlabel='', ylabel="UMI counts")
            g = ax.get_figure()
            g.subplots_adjust(left=0.2)
            g.set_figheight(4)
            g.set_figwidth(4)
            g.savefig(os.path.join(args.output_path, "pathogen_copy_"+virusname_rename+".png"), dpi=500)
            g.savefig(os.path.join(args.output_path, "pathogen_copy_"+virusname_rename+".svg"), dpi=500)

        else:
            print ("STATUS: "+virus_name+" No cells have viral reads not generating violin plot")

def generate_viral_gene_plots(args, dfumidict):
    print ("STATUS: Generating violin plots for viral copies")
    for virus_name, dfumi in dfumidict.items():
        if dfumi is not False:
            sns.set_theme(style="whitegrid")
            # dfumi[gene_name]=[gene_name]*dfumi.shape[0]
            virusname_rename = re.sub("\/", "_", virus_name )
            virusname_rename = re.sub(" ", "_",virusname_rename)
            ax = sns.violinplot(x="gene", y="umi",data=dfumi)
            xlabels = ax.get_xticklabels()
            ax.set_xticklabels(labels=xlabels, rotation=90)
            ax.set(xlabel='genes', ylabel="UMI counts")
            g = ax.get_figure()
            g.subplots_adjust(left=0.2,bottom=0.5)
            g.set_figheight(4)
            g.set_figwidth(6)
            g.savefig(os.path.join(args.output_path, "pathogen_genes_"+virusname_rename+".png"), dpi=500)
            g.savefig(os.path.join(args.output_path, "pathogen_genes_"+virusname_rename+".svg"), dpi=500)

        else:
            print ("STATUS: No cells have viral "+virus_name+" gene read counts not generating violin plot for viral genes")

def generate_coverage_maps(args, dfumidict):
    
    print("STATUS: generating coverage maps")
    files_genome= mr._get_genome_file(args)
    #Get genome length info
    seq = vc.read_FASTA(files_genome[0])

    for virus_name, dfumi in dfumidict.items():
        if dfumi.shape[0] > 0:
            length_genome = len(seq[virus_name])
            krename = re.sub("\/", "_",virus_name )
            krename = re.sub(" ", "_",krename)
            pysam.view("-@", str(args.processors), "-b",  os.path.join(args.output_path, "pathogen_al_mapped_sort.bam"), 
                    virus_name, "-o",  os.path.join(args.output_path, "pathogen_al_mapped_sort_"+krename+".bam"), catch_stdout=False)
            pysam.index(os.path.join(args.output_path, "pathogen_al_mapped_sort_"+krename+".bam"))

            arg=["bamCoverage", "-b", os.path.join(args.output_path, "pathogen_al_mapped_sort_"+krename+".bam"), 
                    "-o", os.path.join(args.output_path, "coverage_"+krename+".bw")]
            ef._run_subprocesses(arg, "STATUS: Generating bigWig file...", "generating bigWig")

            files_gtf = glob.glob(os.path.join(args.path2genome, "*.gtf"))
            if len(files_gtf) ==1:
                print("STATUS: processing gtf for conversion to bed")
                with open(files_gtf[0], 'r') as fin:
                    with open(os.path.join(args.output_path, "genes_"+virus_name+".gtf"), "w") as fout:
                        for line in fin: 
                            larray = line.split("\t")
                            if(larray[2]=="gene" and larray[0]==virus_name):
                                fout.write(line)

                cmd = "gtf2bed < "+os.path.join(args.output_path, "genes_"+virus_name+".gtf")+" > "+os.path.join(args.output_path, "genes_"+virus_name+".bed") 
                r = subprocess.call(cmd, shell=True)  
                ef._check_subprocess_run(r,r, "gtf 2 bed conversion")
                print("STATUS: generating coverage map...")
                argu=["make_tracks_file", "--trackFile",  os.path.join(args.output_path, "genes_"+krename+".bed"), 
                    os.path.join(args.output_path,"coverage_"+krename+".bw"),
                    "-o", os.path.join(args.output_path, "tracks_"+krename+".ini")]
                ef._run_subprocesses(argu, "STATUS: Generating track ini file...", "generating tracks_"+krename+".ini")
                
                process_ini_file(args, os.path.join(args.output_path, "tracks_"+krename+".ini"))

                argu=["pyGenomeTracks", "--tracks", os.path.join(args.output_path,"tracks_"+krename+"edit.ini"),
                    "--region", str(virus_name)+":"+str(1)+"-"+str(length_genome), "--dpi", "200", "--outFileName",
                    os.path.join(args.output_path, krename+"_coveragemap.pdf")]
                ef._run_subprocesses(argu, "STATUS: Generating pdf coverage map...", "generating pdf coverage map")

                argu=["pyGenomeTracks", "--tracks", os.path.join(args.output_path, "tracks_"+krename+"edit.ini"),
                    "--region", str(virus_name)+":"+str(1)+"-"+str(length_genome), "--dpi", "200", "--outFileName",
                    os.path.join(args.output_path, krename+"_coveragemap.png")]
                ef._run_subprocesses(argu, "STATUS: Generating png coverage map...", "generating png coverage map")
            
            else: 
                print("STATUS: number gtf files is "+str(len(files_gtf))+" not incoroporating track info from gtf")

                print("STATUS: generating coverage map...")
                argu=["make_tracks_file", "--trackFile",  os.path.join(args.output_path, "coverage_"+krename+".bw"),
                    "-o", os.path.join(args.output_path, "tracks_"+krename+".ini")]
                ef._run_subprocesses(argu, "STATUS: Generating track ini file...", "generating tracks_"+krename+".ini")
                
                process_ini_file(args, os.path.join(args.output_path, "tracks_"+krename+".ini"))
                
                argu=["pyGenomeTracks", "--tracks", os.path.join(args.output_path, "tracks_"+krename+"edit.ini"),
                    "--region", str(virus_name)+":"+str(1)+"-"+str(length_genome), "--dpi", "200", "--outFileName",
                    os.path.join(args.output_path, krename+"_coveragemap.pdf")]
                ef._run_subprocesses(argu, "STATUS: Generating pdf coverage map...", "generating pdf coverage map")

                argu=["pyGenomeTracks", "--tracks", os.path.join(args.output_path, "tracks_"+krename+"edit.ini"),
                    "--region", str(virus_name)+":"+str(1)+"-"+str(length_genome), "--dpi", "200","--outFileName",
                    os.path.join(args.output_path, krename+"_coveragemap.png")]
                ef._run_subprocesses(argu, "STATUS: Generating png coverage map...", "generating png coverage map")
        else:
            print ("STATUS: "+virus_name+" No cells have pathogen reads reads not generating coverage plot")