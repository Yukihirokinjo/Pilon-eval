#!/bin/bash
#
# Pilon-eval.bash
#
usage_exit() {
	echo "Usage: Pilon-eval.bash  [-i input_contigs/scaffolds(.fasta)] [-1 reads(F)] [-2 reads(R)] [-o output_dir]" 1>&2
	exit 1
}

version=0.8
Pilon_path=`which pilon-1.22.jar`

##Get options

while [ "$#" -gt 0 ]
do
	case "$1" in
		'-v' | '--version' )
			echo "Ver._$version" 
			exit 1
			;;
		'-h' | '--help' )
			usage_exit
			;;
		'-i')
			if  [  -e "$2"  ]; then
				scaffolds="$2" 
				shift 2
			else
				echo "[Error] The input contigs/scaffolds file is not found" 1>&2
				exit 1	 
			fi
			;;
		'-1')
			if  [  -e "$2"  ]; then
				read_F="$2" 
				shift 2
			else
				echo "[Error] The input file is not found" 1>&2
				exit 1	 
			fi
			;;
		'-2')
			if  [  -e "$2"  ]; then
				read_R="$2" 
				shift 2
			else
				echo "[Error] The input file is not found" 1>&2
				exit 1	 
			fi
			;;
		'-o')
			if  [ ! -e "$2"  ]; then
				out_dir="$2" 
				shift 2
			else
				echo "[Error] The output directory is already exist" 1>&2	 
				exit 1	 
			fi
			;;
		'-c')
			if [ -z "$2" ]; then
				echo "PROGRAM: option requires an argument $1" 1>&2
				exit 1
			else
				if  [ `expr "$2" : "[0-9]*$"` -gt 0  ]; then
					Cpu="$2" 
					shift 2
				else
					echo " Argument with option $1 should be an integer " 1>&2
					exit 1
				fi
			fi
			;;
		'-k')
			if [ -z "$2" ]; then
				echo "PROGRAM: option requires an argument $1" 1>&2
				exit 1
			else
				if  [ `expr "$2" : "[0-9]*$"` -gt 0  ]; then
					Kmer="$2" 
					shift 2
				else
					echo " Argument with option $1 should be an integer " 1>&2
					exit 1
				fi
			fi
			;;
		*)
		echo "Invalid option "$1" " 1>&2 
		usage_exit
		;;
	esac
done

if [ -z "$scaffolds" ]; then
	echo "Input file is not specified" 1>&2
	usage_exit
fi
if [ -z "$read_F" ]; then
	echo "Read file is not specified" 1>&2
	usage_exit
fi
if [ -z "$read_R" ]; then
	echo "Read file is not specified" 1>&2
	usage_exit
fi

[ -z "$out_dir" ] && out_dir="Pilon_eval_`date +%Y%m%d`_`date +%H%M%S`"

# Estimate Insert Size
mkdir ${out_dir}
mkdir ${out_dir}/EstInsSize
cp ${scaffolds} ${out_dir}
bowtie2-build -f ${scaffolds} ${out_dir}/tmp_BT2  > ${out_dir}/BT-build.log
seqtk sample -s100 ${read_F} 0.01 > ${out_dir}/EstInsSize/tmp_subread_F.fastq
seqtk sample -s100 ${read_R} 0.01 > ${out_dir}/EstInsSize/tmp_subread_R.fastq
  bowtie2 -x ${out_dir}/tmp_BT2 -1 ${out_dir}/EstInsSize/tmp_subread_F.fastq  -2 ${out_dir}/EstInsSize/tmp_subread_R.fastq \
        -S ${out_dir}/EstInsSize/subBT2.sam  -5 10 -3 10 \
        -I 100 -X 1000 -p ${Cpu:=1} --score-min L,0,0.1  -5 5 -3 10 -D 5 -R 1 -N 0 -L 30 -i S,0,2.5 \
        --no-mixed  --no-discordant --ignore-quals 2> ${out_dir}/EstInsSize/BT2map.log

rm ${out_dir}/EstInsSize/tmp_subread_*.fastq

grep -v "XS:" ${out_dir}/EstInsSize/subBT2.sam > ${out_dir}/EstInsSize/subBT2.unique.sam
samtools stats ${out_dir}/EstInsSize/subBT2.unique.sam > ${out_dir}/EstInsSize/SamStat.txt
grep ^SN ${out_dir}/EstInsSize/SamStat.txt  | cut -f 2-  > ${out_dir}/EstInsSize/SamStatSN.txt

# Set library information (insert size, read length, etc.) 
InsSizeT=`grep "insert size average" ${out_dir}/EstInsSize/SamStatSN.txt | cut -f 2 `
InsSize=${InsSizeT%.*}
InsDevT=`grep "insert size standard deviation" ${out_dir}/EstInsSize/SamStatSN.txt | cut -f 2`
InsDev=${InsDevT%.*}
ReadLen=`grep "average length:" ${out_dir}/EstInsSize/SamStatSN.txt | cut -f 2`
minInsSize=$(( InsSize - 3*InsDev - 10 ))
maxInsSize=$(( InsSize + 3*InsDev + 10 ))
if [ $minInsSize -lt $ReadLen  ] ; then
  minInsSize=$ReadLen
fi
if [ $maxInsSize -gt 1000  ] ; then
  maxInsSize=1000
fi

echo " Min Insert Size = $minInsSize"
echo " Max Insert Size = $maxInsSize"

# Align reads by bowtie2
bowtie2 -x ${out_dir}/tmp_BT2  -1 ${read_F} -2 ${read_R}   -S ${out_dir}/tmp.sam \
        -p ${Cpu:=1} -I ${minInsSize} -X ${maxInsSize} \
	-5 10 -3 10 --score-min L,0,-0.1  -D 5 -R 1 -N 0 -L 30 -i S,0,2.5 \
         > ${out_dir}/BT-map.log

cd ${out_dir}
  samtools view -@ ${Cpu:=1} -bS tmp.sam > tmp.bam
  samtools sort -@ ${Cpu:=1} tmp.bam -o tmp.sort.bam
  samtools index tmp.sort.bam
  rm tmp.sam tmp.bam
cd ..

java -Xmx32G  -jar ${Pilon_path} --genome ${scaffolds} --frags ${out_dir}/tmp.sort.bam \
     --K ${Kmer:=47} --changes --mingap 1 --outdir ${out_dir} --fix all,breaks,amb,circles  --threads ${Cpu:=1}

#END OF FILE
