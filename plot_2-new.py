import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import uuid
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
from collections import Counter
from Bio.Blast import NCBIXML
import matplotlib as mpl
mpl.use('Agg')

#plt.style.use('ggplot')

def run_blast(query_seq, ref_seq):
    unique_id = uuid.uuid4()
    temp_dir = "/tmp"

    query_filename = os.path.join(temp_dir, f"temp_query_{unique_id}.fasta")
    ref_filename = os.path.join(temp_dir, f"temp_ref_{unique_id}.fasta")
    output_filename = os.path.join(temp_dir, f"temp_blast_{unique_id}.xml")

    with open(query_filename, "w") as query_file:
        query_file.write(">query_sequence\n" + query_seq)

    with open(ref_filename, "w") as ref_file:
        ref_file.write(">ref_sequence\n" + ref_seq)

    blastn_cline = NcbiblastnCommandline(query=query_filename, subject=ref_filename, outfmt=5, out=output_filename)
    blastn_cline()

    hsp_intervals = []

    with open(output_filename) as result_handle:
        blast_record = NCBIXML.read(result_handle)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                hsp_intervals.append((hsp.sbjct_start, hsp.sbjct_end))

    os.remove(query_filename)
    os.remove(ref_filename)
    os.remove(output_filename)

    # 合并重叠的 HSP 区间并计算总长度
    if not hsp_intervals:
        return 0  # 或者 return None，看你如何处理无匹配的情况

    hsp_intervals = sorted(hsp_intervals, key=lambda x: x[0])
    merged_intervals = []
    current_start, current_end = hsp_intervals[0]

    for start, end in hsp_intervals[1:]:
        if start <= current_end:  # 如果区间重叠或相邻
            current_end = max(current_end, end)
        else:
            merged_intervals.append((current_start, current_end))
            current_start, current_end = start, end

    merged_intervals.append((current_start, current_end))
    total_coverage_length = sum(end - start + 1 for start, end in merged_intervals)

    return total_coverage_length

def calculate_full_coverage_lengths(input_fastq, ref_sequence):
    new_lengths = {}
    sequence_dict = {record.id: str(record.seq) for record in SeqIO.parse(input_fastq, "fastq")}

    with Pool(processes=cpu_count()) as pool:
        coverage_lengths = pool.starmap(run_blast, [(seq, ref_sequence) for seq in sequence_dict.values()])

    for record_id, coverage_length in zip(sequence_dict.keys(), coverage_lengths):
        new_lengths[record_id] = coverage_length

    return new_lengths

def plot_length_distribution(lengths_dict, output_filename, fasta_len=None, max_length=None):
    if not lengths_dict:
        return

    lengths = list(lengths_dict.values())
    length_counts = Counter(lengths)
    total_sequences = sum(length_counts.values())

    sorted_lengths = sorted(length_counts.keys())
    frequencies = [length_counts[length] for length in sorted_lengths]

    cumulative_sequences = 0
    cumulative_percentages = [0]
    length_percentages = [0]

    for length in sorted_lengths:
        cumulative_sequences += length_counts[length]
        cumulative_percentage = (cumulative_sequences / total_sequences) * 100
        length_percentages.append(length)
        cumulative_percentages.append(cumulative_percentage)
        length_percentages.append(length)
        cumulative_percentages.append(cumulative_percentage + (0.1 / total_sequences))

    length_percentages.append(length_percentages[-1])
    cumulative_percentages.append(100)

    plt.figure(figsize=(10, 6), dpi=300)
    plt.fill_between(cumulative_percentages, length_percentages, step="pre", color='skyblue', alpha=0.7)
    plt.xlabel("Cumulative Percentage (%)", fontsize=20, weight='bold')
    plt.ylabel("Length (nt)", fontsize=20, weight='bold')
    plt.title("Sequence Length Distribution", fontsize=20, weight='bold')
    if max_length:
        plt.ylim(0, max_length * 1.05)

    if fasta_len:
        plt.axhline(y=fasta_len, color='r', lw=0.7, linestyle='--')
        plt.text(0, fasta_len, f'ref_len = {fasta_len}', color='r',
                verticalalignment='bottom', fontsize=12, weight='bold')
    plt.xlim(0, 100)
    plt.grid(False)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.savefig(output_filename)
    plt.close()

def process_files(demultiplexed_folder, reference_folder, output_folder1, output_folder2):
    # 遍历 Demultiplexed 文件夹中的所有 fastq 文件
    for fastq_file in os.listdir(demultiplexed_folder):
        if fastq_file.endswith(".fastq"):
            fastq_path = os.path.join(demultiplexed_folder, fastq_file)
            fasta_file = fastq_file.replace(".fastq", ".fasta")
            ref_path = os.path.join(reference_folder, fasta_file)

            # 检查对应的参考文件是否存在
            if os.path.exists(ref_path):
                print(f"Processing {fastq_file} with reference {fasta_file}...")
                ref_seq = str(SeqIO.read(ref_path, "fasta").seq)
                original_lengths = {record.id: len(record.seq) for record in SeqIO.parse(fastq_path, "fastq")}
                new_lengths = calculate_full_coverage_lengths(fastq_path, ref_seq)
                max_length_ori = max(original_lengths.values()) * 1.5
                original_output_image = os.path.join(output_folder2, fastq_file.replace(".fastq", "_original_LD.png"))
                plot_length_distribution(original_lengths, original_output_image, len(ref_seq), max_length_ori)

                if new_lengths:
                    max_length_new = max(new_lengths.values()) * 1.5
                    output_image = os.path.join(output_folder1, fastq_file.replace(".fastq", "_LD.png"))
                    plot_length_distribution(new_lengths, output_image, len(ref_seq), max_length_new)
                else:
                    print(f"No sequences to plot for {fastq_file}.")
            else:
                print(f"Reference file {fasta_file} not found for {fastq_file}.")

if __name__ == "__main__":
    demultiplexed_folder = "Demultiplexed"  # 输入的FASTQ文件所在文件夹路径
    reference_folder = "Reference"  # 参考序列的FASTA文件所在文件夹路径
    output_folder1 = "净重"  # 保存经过 BLAST 比对后的图像
    output_folder2 = "毛重"  # 保存原始数据的图像

    os.makedirs(output_folder1, exist_ok=True)  # 如果输出文件夹不存在，则创建
    os.makedirs(output_folder2, exist_ok=True)

    process_files(demultiplexed_folder, reference_folder, output_folder1, output_folder2)
