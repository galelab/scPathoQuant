#!/bin/bash

#PBS -N pathoquant
#PBS -l nodes=1:ppn=50
#PBS -l walltime=120:00:00
#PBS -l mem=130gb

# Startup scripts

source $HOME/anaconda3/bin/activate
conda activate scpathoquant

# environment variables

array=("P63T-E" \
"P61T-E" \
"P57T-E" \
"P56T-E" \
"P62T-E" \
"P54T-E" \
"P52T-E" \
"P49T-E" \
"P48T-E" \
"P17T-I" \
"P44T-E" \
"P42T-E" \
"P47T-E" \
"P40T-E" \
"P39T-E" \
"P38T-E" \
"P37T-E" \
"P36T-E" \
"P32T-E" \
"P31T-E" \
"P30T-E" \
"P16T-I" \
"P27T-E" \
"P26T-E" \
"P24T-E" \
"P28T-E" \
"P23T-E" \
"P22T-E" \
"P21T-E" \
"P20T-E" \
"P17T-E" \
"P15T-I" \
"P16T-E" \
"P19T-E" \
"P15T-E" \
"P12T-E" \
"P11T-E" \
"P10T-E" \
"P8T-E" \
"P2T-E" \
"P1T-E" \
"P9T-E" \
"P130N-I" \
"P12T-I" \
"P128N-I" \
"P127N-I" \
"P126N-I" \
"P130T-I" \
"P128T-I" \
"P127T-I" \
"P126T-I" \
"P107T-I" \
"P104T-I" \
"P94T-I" \
"P11T-I" \
"P91T-I" \
"P89T-I" \
"P87T-I" \
"P84T-I" \
"P83T-I" \
"P82T-I" \
"P79T-I" \
"P75T-I" \
"P74T-I" \
"P80T-I" \
"P10T-I" \
"P65T-I" \
"P63T-I" \
"P62T-I" \
"P57T-I" \
"P56T-I" \
"P54T-I" \
"P61T-I" \
"P52T-I" \
"P49T-I" \
"P48T-I" \
"P9T-I" \
"P44T-I" \
"P42T-I" \
"P40T-I" \
"P47T-I" \
"P39T-I" \
"P38T-I" \
"P37T-I" \
"P36T-I" \
"P32T-I" \
"P31T-I" \
"P8T-I" \
"P30T-I" \
"P28T-I" \
"P27T-I" \
"P26T-I" \
"P24T-I" \
"P23T-I" \
"P22T-I" \
"P76-I" \
"P76-E" \
"P130N-E" \
"P128N-E" \
"P21T-I" \
"P126N-E" \
"P130T-E" \
"P128T-E" \
"P127N-E" \
"P127T-E" \
"P126T-E" \
"P107T-E" \
"P104T-E" \
"P91T-E" \
"P20T-I" \
"P89T-E" \
"P94T-E" \
"P87T-E" \
"P84T-E" \
"P83T-E" \
"P82T-E" \
"P79T-E" \
"P75T-E" \
"P74T-E" \
"P80T-E" \
"P65T-E" \
"P19T-I" \
"P2T-I" \
"P1T-I")
array2=("SRR15093488" "SRR15093489" "SRR15093490" "SRR15093491" "SRR15093492"  "SRR15093493" "SRR15093494" "SRR15093495" "SRR15093496" "SRR15093497" "SRR15093498" "SRR15093499" "SRR15093500" "SRR15093501" "SRR15093502" "SRR15093503" "SRR15093504" "SRR15093505" "SRR15093506" "SRR15093507" "SRR15093508" "SRR15093509" "SRR15093510" "SRR15093511" "SRR15093512" "SRR15093513" "SRR15093514" "SRR15093515" "SRR15093516" "SRR15093517" "SRR15093518" "SRR15093519" "SRR15093520" "SRR15093521" "SRR15093522" "SRR15093523" "SRR15093524" "SRR15093525" "SRR15093526" "SRR15093527" "SRR15093528" "SRR15093529" "SRR15093530" "SRR15093531" "SRR15093532" "SRR15093533" "SRR15093534" "SRR15093535" "SRR15093536" "SRR15093537" "SRR15093538" "SRR15093539" "SRR15093540" "SRR15093541" "SRR15093542" "SRR15093543" "SRR15093544" "SRR15093545" "SRR15093546" "SRR15093547" "SRR15093548" "SRR15093549" "SRR15093550" "SRR15093551" "SRR15093552" "SRR15093553" "SRR15093554" "SRR15093555" "SRR15093556" "SRR15093557" "SRR15093558" "SRR15093559" "SRR15093560" "SRR15093561" "SRR15093562" "SRR15093563" "SRR15093564" "SRR15093565" "SRR15093566" "SRR15093567" "SRR15093568" "SRR15093569" "SRR15093570" "SRR15093571" "SRR15093572" "SRR15093573" "SRR15093574" "SRR15093575" "SRR15093576" "SRR15093577" "SRR15093578" "SRR15093579" "SRR15093580" "SRR15093581" "SRR15093582" "SRR15093583" "SRR15093584" "SRR15093585" "SRR15093586" "SRR15093587" "SRR15093588" "SRR15093589" "SRR15093590" "SRR15093591" "SRR15093592" "SRR15093593" "SRR15093594" "SRR15093595" "SRR15093596" "SRR15093597" "SRR15093598" "SRR15093599" "SRR15093600" "SRR15093601" "SRR15093602" "SRR15093603" "SRR15093604" "SRR15093605" "SRR15093606" "SRR15093607" "SRR15093608" "SRR15093609" "SRR15093610" "SRR15093611")
cancer=("Esophagus_cancer")

# Variables
_start=1

# This accounts as the "totalState" variable for the ProgressBar function
_end=124

# Functions

# 1. Create ProgressBar function
# 1.1 Input is currentState($1) and totalState($2)
function ProgressBar {
# Process data
    let _progress=(${1}*100/${2}*100)/100
    let _done=(${_progress}*4)/10
    let _left=40-$_done
# Build progressbar string lengths
    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")

# 1.2 Build progressbar strings and print the ProgressBar line
# 1.2.1 Output example:                           
# 1.2.1.1 Progress : [########################################] 100%
printf "\rProgress : [${_fill// /#}${_empty// /-}] ${_progress}%%"

}

# Script

mkdir -p $HOME/scPathoQuant/results_$cancer/pdf_bucket
for i in "${!array[@]}"; do
	for number in $(seq ${_start} ${_end}); do
		printf "Processing sample %s with SRR code: %s\n" "${array[i]}" "${array2[i]}"
		mkdir -p $HOME/scPathoQuant/results_$cancer/"${array[i]}"/"${array2[i]}"
		scpathoquant -10x $HOME/SC/fastq/$cancer/"${array[i]}"/"${array2[i]}"_S1_L001_ -op $HOME/scPathoQuant/results_$cancer/"${array[i]}"/"${array2[i]}" -p 50 -p2genome $HOME/scPathoQuant/ref_EC/ --tmp_removal
        	cd $HOME/scPathoQuant/results_$cancer/"${array[i]}"/"${array2[i]}"
		PDF_FILE=`ls *.pdf`
		for f in $PDF_FILE; do
			cp $f $HOME/scPathoQuant/results_$cancer/pdf_bucket/"${array2[i]}"_$f
		done
	ProgressBar ${number} ${_end}
	done
done

printf '\nFinished!\n'

exit
