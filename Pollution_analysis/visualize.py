import os
import subprocess
import argparse


def visualizing(reference, fastq_file, outpath):
    """Generate .sorted.bam and .sorted.bam.bai files for a given fastq file."""
    mmifile = f'{os.path.splitext(reference)[0]}.mmi'
    subprocess.run(f'minimap2 -x map-ont -d "{mmifile}" "{reference}"', shell=True)

    purename = os.path.basename(fastq_file).replace(".fastq", "")
    outsam = os.path.join(outpath, f'{purename}.sam')
    outbam = os.path.join(outpath, f'{purename}.bam')
    outsorted_bam = os.path.join(outpath, f'{purename}.sorted.bam')

    # Align, convert to BAM, sort, and index
    subprocess.run(f'minimap2 -ax map-ont "{mmifile}" "{fastq_file}" > "{outsam}"', shell=True)
    subprocess.run(f'samtools view -bS "{outsam}" > "{outbam}"', shell=True)
    subprocess.run(f'samtools sort -O bam -o "{outsorted_bam}" -T temp "{outbam}"', shell=True)
    subprocess.run(f'samtools index "{outsorted_bam}"', shell=True)

    os.remove(outsam)
    os.remove(outbam)
    return outsorted_bam

parser = argparse.ArgumentParser(description="可视化")
parser.add_argument('--fq',  type=str, required=True, help="输入的 FASTQ 文件路径")
parser.add_argument('--ref', type=str, required=True, help="参考序列的 fasta 文件")
parser.add_argument('--out', type=str, required=True, help="输出目录")
args = parser.parse_args()
# Assuming '9mm3.fastq' and the reference file are in the current directory
fastq_file = args.fq
reference = args.ref  # Replace with your actual reference file name
outpath = args.out  # Replace with your desired output directory

# Ensure the output directory exists
os.makedirs(outpath, exist_ok=True)

# Call the visualization function
visualized_file = visualizing(reference, fastq_file, outpath)
