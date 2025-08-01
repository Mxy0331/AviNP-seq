#!/usr/bin/env python3

import os
import subprocess
import argparse
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import yaml  # 新增导入用于解析 YAML 文件
import re
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
from collections import Counter
from Bio.Blast import NCBIXML
import matplotlib as mpl
mpl.use('Agg')
def check_blast_installed():
    """检查BLAST是否安装"""
    try:
        subprocess.run(['blastn', '-version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run(['makeblastdb', '-version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print("Error: BLAST is not installed or not found in PATH.")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: BLAST is not installed or not found in PATH.")
        sys.exit(1)


def create_blast_db(fasta_file, db_name):
    """为给定的FASTA文件创建BLAST数据库"""
    if not os.path.exists(fasta_file):
        print(f"Error: Reference file {fasta_file} does not exist.")
        sys.exit(1)
    # 检查是否已创建数据库
    required_files = [f"{db_name}.{ext}" for ext in ["nhr", "nin", "nsq"]]
    if all(os.path.exists(f) for f in required_files):
        print(f"BLAST数据库已经存在：{db_name}")
        return
    cmd = [
        'makeblastdb',
        '-in', fasta_file,
        '-dbtype', 'nucl',
        '-out', db_name
    ]
    print(f"正在创建BLAST数据库：{fasta_file} 作为 {db_name}...")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error creating BLAST数据库 for {fasta_file}: {e}")
        sys.exit(1)
    print("数据库创建完成。")


def fastq_to_fasta(fastq_file, fasta_file):
    """将FastQ文件转换为FASTA格式"""
    if not os.path.exists(fastq_file):
        print(f"Error: FastQ file {fastq_file} does not exist.")
        sys.exit(1)
    print(f"正在将 {fastq_file} 转换为 {fasta_file}...")
    with open(fastq_file, 'r') as fq, open(fasta_file, 'w') as fa:
        for record in SeqIO.parse(fq, "fastq"):
            SeqIO.write(record, fa, "fasta")
    print("转换完成。")


def run_blast(query_fasta, db_name, output_file):
    """运行BLASTN比对"""
    cmd = [
        'blastn',
        '-query', query_fasta,
        '-db', db_name,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid qstart qend sstart send'
    ]
    print(f"运行BLASTN: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running BLASTN: {e}")
        sys.exit(1)
    print(f"BLASTN结果已保存到 {output_file}。")

def reverse_complement(seq):
    """返回序列的反向互补"""
    return str(Seq(seq).reverse_complement())

def create_reference_sequences(itr_file, vector_file, output_dir):
    # 读取ITR序列
    itr_sequences = {}
    with open(itr_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            itr_sequences[record.id] = str(record.seq)

    # 读取vector序列
    with open(vector_file, 'r') as f:
        vector_sequence = str(SeqIO.read(f, 'fasta').seq)

    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 定义参考序列组合的名称
    combinations = [
        ("Left-flip", "forward", "Right-flip"),
        ("Left-flop", "forward", "Right-flip"),
        ("Left-flip", "forward", "Right-flop"),
        ("Left-flop", "forward", "Right-flop"),
        ("Left-flip", "reverse", "Right-flip"),
        ("Left-flop", "reverse", "Right-flip"),
        ("Left-flip", "reverse", "Right-flop"),
        ("Left-flop", "reverse", "Right-flop")
    ]

    # 处理每个组合，生成参考序列并写入文件
    for left_itr, direction, right_itr in combinations:
        # 获取ITR序列
        left_seq = itr_sequences.get(left_itr)
        right_seq = itr_sequences.get(right_itr)

        # 根据direction决定vector的方向
        if direction == "forward":
            vector_seq = vector_sequence  # 正向vector
        elif direction == "reverse":
            vector_seq = reverse_complement(vector_sequence)  # 反向互补vector
        else:
            continue

        # 拼接参考序列
        reference_sequence = left_seq + vector_seq + right_seq

        # 创建序列记录
        record = SeqRecord(Seq(reference_sequence), id=f"{left_itr}-{direction}-{right_itr}", description="")

        # 写入FASTA文件
        output_file = os.path.join(output_dir, f"{left_itr}_{direction}_{right_itr}.fasta")
        SeqIO.write(record, output_file, "fasta")
        print(f"已保存参考序列: {output_file}")
def parse_blast_output(blast_output, reference_type):
    """
    解析BLAST输出文件，返回一个字典
    key: read_id
    value: list of alignment positions
    """
    alignments = {}
    with open(blast_output, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue  # 跳过格式不正确的行
            qseqid, sseqid, qstart, qend, sstart, send = parts
            if qseqid not in alignments:
                alignments[qseqid] = []
            alignment_info = {
                'reference': reference_type,
                'sseqid': sseqid,
                'sstart': int(sstart),
                'send': int(send),
                'qstart': int(qstart),
                'qend': int(qend)
            }
            alignments[qseqid].append(alignment_info)
    return alignments


def get_reference_lengths(fasta_file):
    """
    读取FASTA文件，返回一个字典
    key: reference_id
    value: reference_length
    """
    ref_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        ref_lengths[record.id] = len(record.seq)
    return ref_lengths


def merge_itr_alignments(itr_alignments, proximity=10):
    """
    对每条序列的ITR比对结果进行合并，
    合并重叠、相连或相隔不超过proximity碱基的区间。
    返回一个字典，key: read_id, value: list of merged [start, end]
    """
    merged = {}
    for read, aln_list in itr_alignments.items():
        # 提取所有的 (qstart, qend) 区间
        intervals = sorted([(aln['qstart'], aln['qend']) for aln in aln_list], key=lambda x: x[0])
        merged_intervals = []
        for interval in intervals:
            if not merged_intervals:
                merged_intervals.append(list(interval))
            else:
                last = merged_intervals[-1]
                # 如果当前区间与上一个区间重叠、相连，或相隔不超过proximity碱基，合并它们
                if interval[0] <= last[1] + proximity:
                    merged_intervals[-1][1] = max(last[1], interval[1])
                else:
                    merged_intervals.append(list(interval))
        merged[read] = merged_intervals
    return merged


def merge_alignments(itr_alignments, vector_alignments, backbone_alignments):
    """合并ITR、Vector和Backbone的比对结果"""
    all_reads = set(list(itr_alignments.keys()) + list(vector_alignments.keys()) + list(backbone_alignments.keys()))
    merged = {}
    for read in all_reads:
        merged[read] = {
            'ITR': itr_alignments.get(read, []),
            'Vector': vector_alignments.get(read, []),
            'Backbone': backbone_alignments.get(read, [])  # 新增Backbone
        }
    return merged


def classify_and_generate_files(itr_flipflop_types, vector_alignments, merged_alignments,  plot_dir,
                                read_lengths, input_file, vector_ref_lengths, backbone_ref_lengths):
    # 定义分类名称
    categories = [
        "flip-minus-flip", "flop-minus-flip", "flip-minus-flop", "flop-minus-flop",
        "flip-plus-flip", "flop-plus-flip", "flop-plus-flop", "flip-plus-flop"
    ]
    qualifying_reads = filter_reads_by_vector_num(
        vector_alignments,
        vector_ref_lengths,
        vec_num=1
    )
    extract_reads(input_file, qualifying_reads, 'single_tmp_file.fastq')
    # 逐分类处理
    for category in categories:
        category_reads = set()

        # 根据分类名称筛选对应的read_ids
        for read_id, alignments in merged_alignments.items():
            itr_types = [itr_flipflop_types.get(read_id, {}).get(0, ''), itr_flipflop_types.get(read_id, {}).get(1, '')]
            vector_orientation = 'forward' if vector_alignments.get(read_id, [{}])[0].get('sstart', 0) <= \
                                              vector_alignments.get(read_id, [{}])[0].get('send', 0) else 'reverse'

            # 判断分类逻辑
            if category == "flip-minus-flip" and itr_types[0] in ['Left-flip', 'Right-flip'] and itr_types[1] in ['Right-flip', 'Left-flip'] and vector_orientation == 'reverse':
                category_reads.add(read_id)
            elif category == "flop-minus-flip" and itr_types[0] in ['Left-flop', 'Right-flop'] and itr_types[1] in ['Right-flip', 'Left-flip'] and vector_orientation == 'reverse':
                category_reads.add(read_id)
            elif category == "flip-minus-flop" and itr_types[0] in ['Left-flip', 'Right-flip'] and itr_types[1] in ['Right-flop', 'Left-flop'] and vector_orientation == 'reverse':
                category_reads.add(read_id)
            elif category == "flop-minus-flop" and itr_types[0] in ['Left-flop', 'Right-flop'] and itr_types[1] in ['Right-flop', 'Left-flop'] and vector_orientation == 'reverse':
                category_reads.add(read_id)
            elif category == "flip-plus-flip" and itr_types[0] in ['Left-flip', 'Right-flip'] and itr_types[1] in ['Right-flip', 'Left-flip'] and vector_orientation == 'forward':
                category_reads.add(read_id)
            elif category == "flop-plus-flip" and itr_types[0] in ['Left-flop', 'Right-flop'] and itr_types[1] in ['Right-flip', 'Left-flip'] and vector_orientation == 'forward':
                category_reads.add(read_id)
            elif category == "flop-plus-flop" and itr_types[0] in ['Left-flop', 'Right-flop'] and itr_types[1] in ['Right-flop', 'Left-flop'] and vector_orientation == 'forward':
                category_reads.add(read_id)
            elif category == "flip-plus-flop" and itr_types[0] in ['Left-flip', 'Right-flip'] and itr_types[1] in ['Right-flop', 'Left-flop'] and vector_orientation == 'forward':
                category_reads.add(read_id)


        # 如果该分类有符合条件的reads
        if category_reads:
            # 计算符合条件的交集
            category_qualifying_reads = set(category_reads) & set(qualifying_reads)

            # 如果交集数量超过100，随机选择100条
            if len(category_qualifying_reads) > 100:
                category_qualifying_reads = random.sample(list(category_qualifying_reads), 100)
                print(f"分类 {category} 的符合条件的序列数量超过100，已随机选择100条进行示意图生成。")

            # 提取FASTQ文件
            fastq_file = f"{category}.fastq"
            extract_reads('single_tmp_file.fastq', category_reads, fastq_file)

            # 生成图像
            plot_all_alignments({read_id: merged_alignments[read_id] for read_id in category_qualifying_reads},
                                read_lengths, plot_dir, vector_ref_lengths, backbone_ref_lengths, itr_flipflop_types,
                                category, title=category)

            print(f"分类 {category} 的 {len(category_qualifying_reads)} 条序列已经生成。")
        else:
            print(f"分类 {category} 没有符合条件的序列。")

def write_sorted_output(merged_alignments, output_file, vector_ref_lengths, backbone_ref_lengths, itr_flipflop_types):
    """
    将合并并排序的比对结果写入输出文件。
    对于每一条序列，按照在测序序列中的比对位置从小到大排序。
    对于ITR和Backbone，只记录测序序列上的位置，并合并重叠区间。
    对于Vector，记录详细的比对信息并计算参考序列上的长度占比。
    """
    with open(output_file, 'w') as out:
        # 写入标题
        out.write("Read_ID\tSorted_Alignments\n")
        for read, alignments in merged_alignments.items():
            alignment_entries = []
            # 处理ITR比对
            itr_merged = alignments['ITR']
            for idx, interval in enumerate(itr_merged):
                # 在此处添加调试信息，检查read和idx是否匹配
                print(f"正在查找 {read} 的第 {idx} 个 ITR类型")

                # 从itr_flipflop_types获取itr_type
                itr_type = itr_flipflop_types.get(read, {}).get(idx, 'Unknown')

                # 如果没有找到itr_type，打印调试信息
                if itr_type == 'Unknown':
                    print(f"未找到 {read} 的第 {idx} 个 ITR类型")

                # 拼接比对信息
                align_str = f"[{interval[0]},{interval[1]}]匹配ITR（{itr_type}）"
                alignment_entries.append({
                    'qstart': interval[0],
                    'entry': align_str
                })

            # 处理Vector比对
            for aln in alignments['Vector']:
                q_range = f"[{aln['qstart']},{aln['qend']}]"
                ref_id = aln['sseqid']
                ref_start = aln['sstart']
                ref_end = aln['send']
                ref_length = vector_ref_lengths.get(ref_id, None)
                if ref_length is None:
                    print(f"Warning: Reference ID {ref_id} not found in Vector reference lengths.")
                    ref_range = f"[{ref_start},{ref_end}]"
                else:
                    ref_start_pct = (ref_start / ref_length) * 100
                    ref_end_pct = (ref_end / ref_length) * 100
                    # 格式化百分比，保留两位小数
                    ref_range = f"[{ref_start} ({ref_start_pct:.2f}%),{ref_end} ({ref_end_pct:.2f}%)]"
                align_str = f"{q_range}匹配Vector {ref_id} {ref_range}"
                alignment_entries.append({
                    'qstart': aln['qstart'],
                    'entry': align_str
                })
            # 处理Backbone比对
            for aln in alignments['Backbone']:
                q_range = f"[{aln['qstart']},{aln['qend']}]"
                ref_id = aln['sseqid']
                ref_start = aln['sstart']
                ref_end = aln['send']
                ref_length = backbone_ref_lengths.get(ref_id, None)
                if ref_length is None:
                    print(f"Warning: Reference ID {ref_id} not found in Backbone reference lengths.")
                    ref_range = f"[{ref_start},{ref_end}]"
                else:
                    ref_start_pct = (ref_start / ref_length) * 100
                    ref_end_pct = (ref_end / ref_length) * 100
                    # 格式化百分比，保留两位小数
                    ref_range = f"[{ref_start} ({ref_start_pct:.2f}%),{ref_end} ({ref_end_pct:.2f}%)]"
                align_str = f"{q_range}匹配Backbone {ref_id} {ref_range}"
                alignment_entries.append({
                    'qstart': aln['qstart'],
                    'entry': align_str
                })
            if not alignment_entries:
                continue  # 如果没有比对信息，跳过
            # 按照qstart排序
            sorted_alignments = sorted(alignment_entries, key=lambda x: x['qstart'])
            # 提取排序后的比对字符串
            sorted_align_strs = [item['entry'] for item in sorted_alignments]
            # 合并为一个字符串，用分号分隔
            merged_align_str = "; ".join(sorted_align_strs)
            # 写入文件
            out.write(f"{read}\t{merged_align_str}\n")
    print(f"已将排序和格式化后的比对结果写入 {output_file}。")


def filter_reads_by_vector_num(vector_alignments, vector_ref_lengths, coverage_threshold=0.9, vec_num=1):
    """
    筛选出包含指定数量（vec_num）完整Vector匹配的序列。
    完整Vector匹配定义为覆盖率≥coverage_threshold。
    返回一个集合，包含符合条件的read_ids。
    """
    qualifying_reads = set()
    for read, aln_list in vector_alignments.items():
        complete_vector_matches = 0
        for aln in aln_list:
            ref_id = aln['sseqid']
            ref_length = vector_ref_lengths.get(ref_id, None)
            if ref_length is None:
                print(f"Warning: Reference ID {ref_id} not found in Vector reference lengths.")
                continue
            # 计算alignment_length为绝对值，确保正向和反向都能正确计算
            alignment_length = abs(aln['send'] - aln['sstart']) + 1
            coverage = alignment_length / ref_length
            if coverage >= coverage_threshold:
                complete_vector_matches += 1
        if complete_vector_matches == vec_num:
            qualifying_reads.add(read)
    return qualifying_reads

def run_bwa_index(reference_file):
    """生成参考序列的BWA索引"""
    if not os.path.exists(reference_file):
        print(f"Error: {reference_file} does not exist.")
        return
    print(f"正在为 {reference_file} 生成BWA索引...")
    cmd = ['bwa', 'index', reference_file]
    subprocess.run(cmd, check=True)
    print(f"参考序列 {reference_file} 的BWA索引生成完成。")

def run_bwa_mem(fastq_file, reference_file, output_sam):
    """使用BWA进行比对，生成SAM文件"""
    if not os.path.exists(fastq_file):
        print(f"Error: FastQ file {fastq_file} does not exist.")
        return
    print(f"正在使用BWA进行比对：{fastq_file} 与 {reference_file}...")
    cmd = ['bwa', 'mem', reference_file, fastq_file]
    with open(output_sam, 'w') as sam_file:
        subprocess.run(cmd, check=True, stdout=sam_file)
    print(f"BWA比对完成，结果已保存到 {output_sam}。")

def sam_to_bam(sam_file, bam_file):
    """将SAM文件转换为BAM文件"""
    print(f"正在将 {sam_file} 转换为BAM文件...")
    cmd = ['samtools', 'view', '-bS', sam_file, '-o', bam_file]
    subprocess.run(cmd, check=True)
    print(f"BAM文件已保存到 {bam_file}。")

def sort_bam(bam_file, sorted_bam_file):
    """对BAM文件进行排序"""
    print(f"正在排序BAM文件 {bam_file}...")
    cmd = ['samtools', 'sort', bam_file, '-o', sorted_bam_file]
    subprocess.run(cmd, check=True)
    print(f"BAM文件已排序并保存为 {sorted_bam_file}。")

def index_bam(sorted_bam_file):
    """为排序后的BAM文件生成索引"""
    print(f"正在为排序后的BAM文件 {sorted_bam_file} 生成索引...")
    cmd = ['samtools', 'index', sorted_bam_file]
    subprocess.run(cmd, check=True)
    print(f"BAM文件索引已生成：{sorted_bam_file}.bai")

def process_fastq_files(fastq_files, reference_files, output_dir):
    """为每个fastq文件与参考序列生成BAM文件"""
    os.makedirs(output_dir, exist_ok=True)

    # 生成参考序列索引文件
    for reference_file in reference_files:
        run_bwa_index(reference_file)

    # 处理每个fastq文件并生成BAM文件
    for fastq_file, reference_file in zip(fastq_files, reference_files):
        if not os.path.exists(fastq_file):
            print(f"警告: FastQ 文件 {fastq_file} 不存在，跳过该文件。")
            continue

        base_name = os.path.splitext(os.path.basename(fastq_file))[0]

        # 构建SAM文件、BAM文件、排序后的BAM文件路径
        output_sam = os.path.join(output_dir, f"{base_name}.sam")
        output_bam = os.path.join(output_dir, f"{base_name}.bam")
        sorted_bam = os.path.join(output_dir, f"{base_name}.sorted.bam")

        # 使用BWA进行比对
        run_bwa_mem(fastq_file, reference_file, output_sam)

        # 将SAM文件转换为BAM文件
        sam_to_bam(output_sam, output_bam)

        # 对BAM文件进行排序
        sort_bam(output_bam, sorted_bam)

        # 为BAM文件生成索引
        index_bam(sorted_bam)

        # 清理中间文件
        os.remove(output_sam)
        os.remove(output_bam)

        print(f"所有步骤完成：{fastq_file} 的BAM文件已生成并排序：{sorted_bam}。\n")

def extract_reads(fastq_file, read_ids, output_fastq):
    """
    从FastQ文件中提取指定的read_ids，并写入到新的FastQ文件中。
    """
    if os.path.abspath(fastq_file) == os.path.abspath(output_fastq):
        print(f"Error: 输出 FastQ 文件 {output_fastq} 与输入 FastQ 文件 {fastq_file} 相同，可能会导致文件被覆盖。")
        sys.exit(1)

    print(f"正在从 {fastq_file} 中提取符合条件的序列到 {output_fastq}...")
    count = 0
    with open(fastq_file, 'r') as infile, open(output_fastq, 'w') as outfile:
        for record in SeqIO.parse(infile, "fastq"):
            if record.id in read_ids:
                SeqIO.write(record, outfile, "fastq")
                count += 1
    print(f"提取完成，共提取了 {count} 条序列到 {output_fastq}。")


def plot_all_alignments(merged_alignments, read_lengths, output_dir, vector_ref_lengths, backbone_ref_lengths,
                        itr_flipflop_types, name, title):
    """
    生成所有序列比对的综合示意图，并保存为PNG文件。
    每个read占据图中的一条水平线。
    ITR匹配部分用黑色表示，位于下方。
    Vector匹配部分用渐变色表示，颜色基于参考序列的位置（红色代表0位，蓝色代表最后一位）。
    Backbone匹配部分用绿色表示，并标识其起始和结束位置。
    同时显示坐标轴和图例。
    """
    fig, ax = plt.subplots(figsize=(20, max(3, len(merged_alignments) * 0.3)))  # 根据read数量调整高度

    # 根据 read_lengths 进行排序，长的在下，短的在上
    sorted_reads = sorted(merged_alignments.keys(), key=lambda read_id: read_lengths.get(read_id, 0))
    current_lengths = [read_lengths[rid] for rid in merged_alignments if rid in read_lengths]
    max_len = max(current_lengths) if current_lengths else 1000

    # 更新y_positions，保证较长的read在下方
    y_positions = {read_id: idx for idx, read_id in enumerate(sorted_reads)}

    ax.set_ylim(-1, len(sorted_reads))
    ax.set_xlim(0, max_len + 1000)  # 动态设置 x 轴上限
    ax.set_xticks(np.arange(0, max_len + 1001, step=1500))
    ax.set_xlabel('Position in Read', fontsize=12)
    ax.set_ylabel('Reads', fontsize=12)
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(sorted_reads, fontsize=6)  # 根据需要调整字体大小
    ax.set_title(f'plot of {title} reads', fontsize=16, loc='center')

    # 添加网格
    ax.grid(axis='x', linestyle='--', alpha=0.5)

    for read_id, alignments in merged_alignments.items():
        y = y_positions[read_id]  # 根据y_positions设置y坐标
        read_length = read_lengths.get(read_id, None)
        if read_length is None:
            print(f"Warning: Read ID {read_id} not found in FastQ file.")
            continue

        ITR_COLORS = {
            'Left-flip': '#FF6B6B',
            'Left-flop': '#4ECDC4',
            'Right-flip': '#f9d580',
            'Right-flop': '#bdb5e1',
            'Unknown': 'gray'
        }

        # 记录所有比对的区间
        all_intervals = []

        for itr_idx, (start, end) in enumerate(alignments['ITR']):
            itr_type = itr_flipflop_types.get(read_id, {}).get(itr_idx, 'Unknown')
            color = ITR_COLORS.get(itr_type, 'gray')
            rect = Rectangle((start, y - 0.2), end - start, 0.4, facecolor=color)
            ax.add_patch(rect)
            all_intervals.append((start, end))

        # 处理Vector比对
        for vector in alignments['Vector']:
            q_start = vector['qstart']
            q_end = vector['qend']
            s_start = vector['sstart']
            s_end = vector['send']
            ref_id = vector['sseqid']
            ref_length = vector_ref_lengths.get(ref_id, None)

            if ref_length is None:
                print(f"Warning: Reference ID {ref_id} not found in Genome reference lengths.")
                continue

            orientation = 'forward' if s_start <= s_end else 'reverse'

            # 计算参考序列的对应位置
            if orientation == 'forward':
                ref_positions = np.linspace(s_start, s_end, 100)
            else:
                ref_positions = np.linspace(s_start, s_end, 100)

            # 标准化参考位置到 [0, 1]
            ref_positions_normalized = (ref_positions - 0) / (ref_length - 1)  # 0-based

            # 生成颜色，根据参考位置
            cmap = cm.get_cmap('coolwarm')  # 使用matplotlib.cm获取colormap
            colors = cmap(ref_positions_normalized)

            # 计算渐变的宽度
            gradient_width = (q_end - q_start) / 100

            for i in range(100):
                rect = Rectangle((q_start + i * gradient_width, y - 0.2), gradient_width, 0.4, facecolor=colors[i],
                                 edgecolor='none')
                ax.add_patch(rect)
            all_intervals.append((q_start, q_end))

        # 处理Backbone比对
        for backbone in alignments['Backbone']:
            q_start = backbone['qstart']
            q_end = backbone['qend']
            s_start = backbone['sstart']
            s_end = backbone['send']
            ref_id = backbone['sseqid']
            ref_length = backbone_ref_lengths.get(ref_id, None)

            if ref_length is None:
                print(f"Warning: Reference ID {ref_id} not found in Backbone reference lengths.")
                continue

            # 使用绿色表示Backbone
            rect = Rectangle((q_start, y - 0.2), q_end - q_start, 0.4, facecolor='green', edgecolor='black')
            ax.add_patch(rect)
            all_intervals.append((q_start, q_end))

        # 标注未匹配段
        all_intervals = sorted(all_intervals, key=lambda x: x[0])
        for i in range(1, len(all_intervals)):
            prev_end = all_intervals[i-1][1]
            curr_start = all_intervals[i][0]
            if curr_start - prev_end > 100:
                rect = Rectangle((prev_end, y - 0.2), curr_start - prev_end, 0.4, facecolor='#FCB8B8', edgecolor='black', hatch='//')
                ax.add_patch(rect)

    # 添加图例
    vector_patch = Rectangle((0, 0), 1, 1, facecolor=cm.get_cmap('coolwarm')(0.0), edgecolor='none')  # 使用cm.get_cmap
    backbone_patch = Rectangle((0, 0), 1, 1, facecolor='green', edgecolor='black')
    ITR_patch_1 = Rectangle((0, 0), 1, 1, facecolor='#FF6B6B', edgecolor='black')
    ITR_patch_2 = Rectangle((0, 0), 1, 1, facecolor='#4ECDC4', edgecolor='black')
    ITR_patch_3 = Rectangle((0, 0), 1, 1, facecolor='#f9d580', edgecolor='black')
    ITR_patch_4 = Rectangle((0, 0), 1, 1, facecolor='#bdb5e1', edgecolor='black')
    unmatched_patch = Rectangle((0, 0), 1, 1, facecolor='#FCB8B8', edgecolor='black', hatch='//')
    ax.legend([vector_patch, backbone_patch, ITR_patch_1, ITR_patch_2, ITR_patch_3, ITR_patch_4, unmatched_patch],
              ['Genome', 'Backbone', 'Left-flip', 'Left-flop', 'Right-flip', 'Right-flop', 'Unmatched'],
              loc='upper right', fontsize=12)

    # 添加颜色条（Colorbar）来表示Vector参考序列的位置
    sm = cm.ScalarMappable(cmap='coolwarm', norm=mcolors.Normalize(vmin=0, vmax=1))  # 使用cm.ScalarMappable
    sm.set_array([])  # 这里传递一个空的数组来初始化颜色条
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('Genome Reference Position', rotation=270, labelpad=15, fontsize=12)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(['Start', 'End'])

    # 保存图像
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, f"{name}.png")
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"示意图已保存到 {output_path}。")

def extract_itr_sequences(fastq_file, merged_itr, output_dir='temp_itr'):
    os.makedirs(output_dir, exist_ok=True)
    temp_files = []
    # read_sequences = SeqIO.to_dict(SeqIO.parse(fastq_file, "fastq"))
    read_sequences = {}
    for record in SeqIO.parse(fastq_file, "fastq"):
        if record.id not in read_sequences:
            read_sequences[record.id] = record

        else:
            # 重复的 read_id 会被忽略并打印警告
            print(f"Warning: Duplicate read ID found and skipped: {record.id}")

    for read_id, intervals in merged_itr.items():
        # 生成安全的文件名前缀（不影响序列ID）
        safe_read_id = re.sub(r'[^a-zA-Z0-9_]', '_', read_id)

        for idx, (start, end) in enumerate(intervals):
            # 生成临时文件名（使用安全前缀）
            temp_file = os.path.join(output_dir, f"{safe_read_id[:20]}_i{idx}.fa")

            if read_id not in read_sequences:
                print(f"Warning: {read_id} not found in FastQ, skipped")
                continue

            seq = str(read_sequences[read_id].seq)
            if end > len(seq):
                print(f"Warning: {read_id} ITR区间[{start},{end}]越界，已跳过")
                continue

            # 关键修复：序列ID使用原始read_id，确保BLAST输出能匹配
            with open(temp_file, 'w') as f:
                f.write(f">{read_id}_itr{idx}\n{seq[start-1:end]}\n")
            temp_files.append(temp_file)

    return temp_files


def run_itr_flipflop_blast(temp_files, db_name, threads, output_file='itr_flipflop_blast.out'):
    """批量运行BLAST并获取构型（Windows兼容版本）"""

    # 打开输出文件，确保每次都以追加模式写入
    with open(output_file, 'a') as f:
        for query_file in temp_files:
            cmd = [
                'blastn',
                '-query', query_file,  # 使用单个查询文件
                '-db', db_name,
                '-outfmt', '6 qseqid sseqid score sstart send',  # 选择输出格式
                '-num_threads', str(threads)
            ]

            # 运行blastn命令，并获取输出结果
            result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # 获取标准输出并逐行处理
            result_output = result.stdout.decode('utf-8')
            for line in result_output.splitlines():
                # 跳过空行
                if line.strip() == '':
                    continue

                # 分割行并检查长度
                parts = line.strip().split('\t')
                if len(parts) != 5:
                    print(f"Warning: Skipping invalid line: {line}")
                    continue

                # 将有效的比对结果写入文件
                f.write(line + "\n")

            print(f"BLAST结果已处理并追加到 {output_file}。")

    print(f"所有BLAST结果已经追加到 {output_file} 中。")



    # 删除临时列表文件
    # os.remove(list_file)

    # 解析结果：{(read_id, itr_index): best_type}
    itr_types = {}
    with open(output_file) as f:
        for line in f:
            qseqid, sseqid, score, sstart, send = line.strip().split('\t')
            if int(sstart) > int(send):
                continue  # 仅保留正向比对
            read_id, itr_idx = qseqid.rsplit('_itr', 1)
            itr_idx = int(itr_idx)
            score = float(score)

            current_key = (read_id, itr_idx)
            if current_key not in itr_types or score > itr_types[current_key]['score']:
                itr_types[current_key] = {'type': sseqid, 'score': score}

    # 转换为最终结构：{read_id: {itr_index: type}}
    final_types = defaultdict(dict)
    for (read_id, itr_idx), data in itr_types.items():
        final_types[read_id][itr_idx] = data['type']
    return final_types


def get_read_lengths(fastq_file):
    """
    读取FastQ文件，返回一个字典
    key: read_id
    value: sequence_length
    """
    read_lengths = {}
    with open(fastq_file, 'r') as fq:
        for record in SeqIO.parse(fq, "fastq"):
            read_lengths[record.id] = len(record.seq)
    return read_lengths


def write_fasta_from_yaml(ref_yaml, output_dir='references'):
    """
    从YAML文件中读取参考序列，并写入单独的FASTA文件。
    返回一个字典，key: reference_type, value: fasta_file_path
    """
    if not os.path.exists(ref_yaml):
        print(f"Error: Reference YAML file {ref_yaml} does not exist.")
        sys.exit(1)

    with open(ref_yaml, 'r') as f:
        try:
            ref_data = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(f"Error parsing YAML file: {exc}")
            sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fasta_files = {}
    for ref_type, ref_info in ref_data.items():
        sequence = ref_info.get('sequence', '').strip().replace('\n', '')
        if not sequence:
            print(f"Warning: No sequence found for {ref_type} in YAML file.")
            continue
        fasta_filename = f"{ref_type}.fasta"
        fasta_path = os.path.join(output_dir, fasta_filename)
        with open(fasta_path, 'w') as fasta_file:
            fasta_file.write(f">{ref_type}\n{sequence}\n")
        fasta_files[ref_type] = fasta_path
        print(f"已从 YAML 文件写入 {ref_type} 到 {fasta_path}。")
    return fasta_files

def get_reads_matched_to_vector(vector_alignments, vector_ref_lengths, coverage_threshold=0.9):
    """
    获取所有完全匹配上vector的reads。
    只筛选那些覆盖率达到给定阈值的reads。

    Args:
        vector_alignments (dict): 含有所有reads与vector的比对信息。
        vector_ref_lengths (dict): vector参考序列的长度字典，key为参考序列ID，value为该序列的长度。
        coverage_threshold (float): 最小的覆盖率阈值，用于筛选。

    Returns:
        set: 包含所有匹配上vector的read_id的集合。
    """
    matched_reads = set()

    for read_id, aln_list in vector_alignments.items():
        complete_vector_matches = 0
        for aln in aln_list:
            ref_id = aln['sseqid']
            ref_length = vector_ref_lengths.get(ref_id, None)

            if ref_length is None:
                print(f"Warning: Reference ID {ref_id} not found in Vector reference lengths.")
                continue

            # 计算alignment_length为绝对值，确保正向和反向都能正确计算
            alignment_length = abs(aln['send'] - aln['sstart']) + 1
            coverage = alignment_length / ref_length

            # 如果覆盖率满足阈值，认为该比对是有效的
            if coverage >= coverage_threshold:
                complete_vector_matches += 1

        # 如果该read符合完全匹配条件，加入结果集
        if complete_vector_matches > 0:
            matched_reads.add(read_id)

    return matched_reads


def plot_length_distribution(qualifying_reads, read_lengths, output_dir, plot_name, fasta_len):
    """
    Plot the cumulative percentage distribution of sequence lengths and save the plot with a specified filename.
    """
    lengths = [read_lengths[read_id] for read_id in qualifying_reads if read_id in read_lengths]
    length_counts = Counter(lengths)
    total_sequences = sum(length_counts.values())  # 总序列数

    sorted_lengths = sorted(length_counts.keys())
    cumulative_sequences = 0
    cumulative_percentages = []
    length_percentages = []

    for length in sorted_lengths:
        cumulative_sequences += length_counts[length]
        cumulative_percentage = (cumulative_sequences / total_sequences) * 100
        cumulative_percentages.append(cumulative_percentage)
        length_percentages.append(length)

    plt.figure(figsize=(10, 6), dpi=300)
    plt.fill_between(cumulative_percentages, length_percentages, step="pre", color='skyblue', alpha=0.8)
    plt.xlabel("Cumulative Percentage (%)")
    plt.ylabel("Length")
    plt.title("Sequence Length Distribution as Cumulative Percentage")
    plt.axhline(y=fasta_len, color='r', lw=0.7, linestyle='--')
    plt.text(0, fasta_len, f'ref_len = {fasta_len}', color='r', verticalalignment='bottom')
    plt.ylim(None, fasta_len * 3)
    plt.grid(True)
    output_path = os.path.join(output_dir, f"{plot_name}.png")
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

    print(f"长度分布图已保存到 {output_path}。")




def parse_arguments():
    parser = argparse.ArgumentParser(
        description="使用BLAST比对FastQ文件到ITR、Vector和Backbone参考序列，并根据指定条件筛选符合条件的序列。")
    # 这里不再需要vec_num和s参数，因为我们将自动生成
    parser.add_argument('-ref', '--ref_yaml', required=True, help='参考序列的 YAML 文件')
    parser.add_argument('-fq', '--fastq_file', required=True, help='输入的FastQ文件')
    parser.add_argument('-o', '--output', default='alignment_results.tsv', help='输出比对结果文件')
    parser.add_argument("-d", "--output_dir", default='plot', help="输出文件的文件夹路径")
    parser.add_argument('-c', '--coverage_threshold', type=float, default=0.9,
                        help='Vector匹配的最小覆盖率阈值（默认0.9，即90%）')
    parser.add_argument('-p', '--proximity', type=int, default=15, help='合并ITR比对区间时的最大间隔（默认10碱基）')
    # 新增以下两个参数
    parser.add_argument('--itr_flipflop', required=True, help='ITR构型参考序列FASTA文件')
    parser.add_argument('--threads', type=int, default=4, help='BLAST比对线程数（默认4）')
    args = parser.parse_args()



    return parser.parse_args()


def generate_statistics(input_file, merged_alignments, vector_alignments, backbone_blast_out='blast_backbone.out',
                        output_dir='statistics'):
    """生成统计Excel文件，包含reads数、总碱基数、各种分类的占比"""

    # 获取文件名作为AAV名称（即FastQ文件的文件名）
    aav_name = os.path.splitext(os.path.basename(input_file))[0]

    # 获取输入FastQ文件的reads数和总碱基数
    total_reads = sum(1 for _ in SeqIO.parse(input_file, "fastq"))
    total_bases = sum(len(record.seq) for record in SeqIO.parse(input_file, "fastq"))

    # 初始化统计数据（所有值均用列表保存，确保长度一致）
    statistics = {
        'AAV Name': [aav_name],
        'Total Reads': [total_reads],
        'Total Bases': [total_bases],
        'Average Read Length': [total_bases / total_reads] if total_reads > 0 else [0]
    }

    # 统计含Backbone的reads数及总对齐碱基数
    backbone_reads = 0
    total_alignment_length = 0  # 存储所有对齐区间的总长度
    with open(backbone_blast_out, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            qstart = int(parts[2])  # 查询序列起始位置
            qend = int(parts[3])    # 查询序列结束位置
            alignment_length = abs(qend - qstart)
            total_alignment_length += alignment_length
            backbone_reads += 1
    backbone_alignment_ratio = total_alignment_length / total_bases if total_bases > 0 else 0

    statistics['Backbone Reads'] = [backbone_reads]
    statistics['Backbone Alignment Bases'] = [total_alignment_length]
    statistics['Backbone Alignment Ratio'] = [backbone_alignment_ratio]

    # 分类统计：0倍、单倍、flip-minus-flip、flop-minus-flip等
    categories = [
        "zero", "single", "flip-minus-flip", "flop-minus-flip", "flip-minus-flop", "flop-minus-flop",
        "flip-plus-flip", "flop-plus-flip", "flop-plus-flop", "flip-plus-flop", "double", "triple", "multi"
    ]

    for category in categories:
        fastq_file_cat = f"{category}.fastq"
        if os.path.exists(fastq_file_cat):
            category_reads = sum(1 for _ in SeqIO.parse(fastq_file_cat, "fastq"))
        else:
            category_reads = 0
        category_ratio = category_reads / total_reads if total_reads > 0 else 0
        statistics[f'{category} Reads'] = [category_reads]
        statistics[f'{category} Ratio'] = [category_ratio]

    # 计算 Complete percentage
    statistics['Complete percentage'] = [((statistics.get('double Reads', [0])[0] +
                                             statistics.get('triple Reads', [0])[0] +
                                             statistics.get('multi Reads', [0])[0] +
                                             statistics.get('single Reads', [0])[0]) / total_reads * 100)
                                          if total_reads > 0 else 0]

    # 计算各类比例时，先取出'single Ratio'
    single_ratio = statistics.get('single Ratio', [0])[0]
    if single_ratio == 0:
        statistics['flip-flip Ratio'] = [0]
        statistics['flip-flop Ratio'] = [0]
        statistics['flop-flop Ratio'] = [0]
        statistics['flop-flip Ratio'] = [0]
    else:
        statistics['flip-flip Ratio'] = [(statistics.get('flip-minus-flip Ratio', [0])[0] +
                                          statistics.get('flip-plus-flip Ratio', [0])[0]) / single_ratio]
        statistics['flip-flop Ratio'] = [(statistics.get('flip-minus-flop Ratio', [0])[0] +
                                          statistics.get('flip-plus-flop Ratio', [0])[0]) / single_ratio]
        statistics['flop-flop Ratio'] = [(statistics.get('flop-minus-flop Ratio', [0])[0] +
                                          statistics.get('flop-plus-flop Ratio', [0])[0]) / single_ratio]
        statistics['flop-flip Ratio'] = [(statistics.get('flop-minus-flip Ratio', [0])[0] +
                                          statistics.get('flop-plus-flip Ratio', [0])[0]) / single_ratio]

    statistics['Total Ratio'] = [statistics['flip-flip Ratio'][0] +
                                 statistics['flip-flop Ratio'][0] +
                                 statistics['flop-flop Ratio'][0] +
                                 statistics['flop-flip Ratio'][0]]

    # 将统计数据转化为DataFrame
    df = pd.DataFrame(statistics)

    # 生成Excel文件
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, f"{aav_name}_statistics.xlsx")
    df.to_excel(output_path, index=False)
    print(f"统计文件已生成并保存为 {output_path}。")


def main():
    args = parse_arguments()

    input_file = args.fastq_file
    output_file = args.output
    plot_dir = args.output_dir


    print(f"开始处理文件：{input_file}")
    vec_nums = [0, 1, 2, 3, 4]
    # 检查BLAST是否安装
    check_blast_installed()

    # 从 YAML 文件中读取参考序列并写入FASTA文件
    references = write_fasta_from_yaml(args.ref_yaml)

    # 确保所有必要的参考序列都已提供
    required_refs = ['ITR', 'vector', 'backbone']
    for ref in required_refs:
        if ref not in references:
            print(f"Error: 参考序列 {ref} 未在 YAML 文件中提供。")
            sys.exit(1)

    # 获取Vector和Backbone参考序列的长度
    ITR_ref_lengths = get_reference_lengths(references['ITR'])
    vector_ref_lengths = get_reference_lengths(references['vector'])
    backbone_ref_lengths = get_reference_lengths(references['backbone'])  # 获取Backbone参考序列长度

    # 创建BLAST数据库
    create_blast_db(references['ITR'], 'ITR_db')
    create_blast_db(references['vector'], 'Vector_db')
    create_blast_db(references['backbone'], 'Backbone_db')  # 创建Backbone BLAST数据库


    # 转换FastQ到Fasta
    fastq_fasta = 'input_sequences.fasta'
    if os.path.abspath(fastq_fasta) == os.path.abspath(input_file):
        print(f"Error: 输出FASTA文件 {fastq_fasta} 与输入FastQ文件 {input_file} 相同，可能会导致文件被覆盖。")
        sys.exit(1)
    fastq_to_fasta(input_file, fastq_fasta)

    # 运行BLASTN比对到ITR、Vector和Backbone
    itr_blast_out = 'blast_itr.out'
    vector_blast_out = 'blast_vector.out'
    backbone_blast_out = 'blast_backbone.out'  # 新增Backbone BLAST输出文件

    run_blast(fastq_fasta, 'ITR_db', itr_blast_out)
    run_blast(fastq_fasta, 'Vector_db', vector_blast_out)
    run_blast(fastq_fasta, 'Backbone_db', backbone_blast_out)  # 运行Backbone BLAST

    # 解析BLAST输出
    itr_alignments_raw = parse_blast_output(itr_blast_out, 'ITR')
    vector_alignments = parse_blast_output(vector_blast_out, 'Vector')
    backbone_alignments = parse_blast_output(backbone_blast_out, 'Backbone')  # 解析Backbone BLAST输出
    # print(f'111111111111{backbone_alignments}')

    # 合并ITR比对区间，设置相邻区间的最大间隔为args.proximity碱基
    itr_alignments = merge_itr_alignments(itr_alignments_raw, proximity=args.proximity)
    total_reads = sum(1 for _ in SeqIO.parse(input_file, "fastq"))
    temp_files = extract_itr_sequences(input_file, itr_alignments)

    if not temp_files:
        print("警告：没有有效的ITR区间需要分析构型")
        itr_flipflop_types = {}
    else:
        # Step 2: 运行ITR构型BLAST
        create_blast_db(args.itr_flipflop, 'ITR_flipflop_db')  # 确保数据库存在
        itr_flipflop_types = run_itr_flipflop_blast(
            temp_files,
            'ITR_flipflop_db',
            args.threads
        )

    # 合并比对结果
    merged_alignments = merge_alignments(itr_alignments, vector_alignments, backbone_alignments)  # 包含Backbone

    # 写入排序和格式化后的输出文件，包含Vector和Backbone参考序列的长度占比
    write_sorted_output(merged_alignments, output_file, vector_ref_lengths, backbone_ref_lengths, itr_flipflop_types)  # 传递Backbone长度
    sample_reads = set()
    for vec_num in vec_nums:
        # 筛选符合条件的read_ids
        if vec_num == 0:
            fastq_name = "zero.fastq"
            plot_name = "zero"
            title='Incomplete Sequence (IS)'
        elif vec_num == 1:
            fastq_name = "single.fastq"
            plot_name = "single"
            title = 'Single Full-Length Sequence (S-FL)'
            read_lengths = get_read_lengths(input_file)
            qualifying_reads = filter_reads_by_vector_num(
                vector_alignments,
                vector_ref_lengths,
                coverage_threshold=args.coverage_threshold,
                vec_num=vec_num

            )
            extract_reads(input_file, qualifying_reads, fastq_name)
            classify_and_generate_files(itr_flipflop_types, vector_alignments, merged_alignments, plot_dir, read_lengths, input_file, vector_ref_lengths, backbone_ref_lengths)
        elif vec_num == 2:
            fastq_name = "double.fastq"
            plot_name = "double"
            title = 'Double Full-Length Sequences (D-FL)'
        elif vec_num == 3:
            fastq_name= "triple.fastq"
            plot_name = "triple"
            title = 'Triple Full-Length Sequences (T-FL)'
        else:
            fastq_name = "multi.fastq"
            plot_name = "multi"
            title = 'Multiple Full-Length Sequences (M-FL)'
        qualifying_reads = filter_reads_by_vector_num(
            vector_alignments,
            vector_ref_lengths,
            coverage_threshold=args.coverage_threshold,
            vec_num=vec_num

        )
        print(f"共识别出 {len(qualifying_reads)} 条包含{vec_num}个完整Vector匹配的序列。")

            # plot_all_alignments(merged_alignments, read_lengths, plot_dir, vector_ref_lengths, backbone_ref_lengths,itr_flipflop_types, f"{plot_name}_classified", categorized_reads)

        # 提取符合条件的序列到新的FastQ文件
        extract_reads(input_file, qualifying_reads, fastq_name)

        # 生成综合示意图


        if qualifying_reads:
            print("正在生成示意图...")
            sample_size = min(int(len(qualifying_reads) / total_reads * 100),len(qualifying_reads))
            sample_reads.update(random.sample(list(qualifying_reads), sample_size))
            # 如果qualifying_reads数量大于100，随机选择100条
            if len(qualifying_reads) > 100:
                qualifying_reads = random.sample(list(qualifying_reads), 100)
                print("条数超过100，已随机选择100条进行示意图生成。")

            # 获取所有read的长度
            read_lengths = get_read_lengths(input_file)

            # 准备所有符合条件的reads的数据
            combined_alignments = {read_id: merged_alignments[read_id] for read_id in qualifying_reads if
                                   read_id in merged_alignments}

            # 生成示意图
            plot_all_alignments(combined_alignments, read_lengths, plot_dir, vector_ref_lengths, backbone_ref_lengths,
                                itr_flipflop_types, plot_name, title)  # 传递Backbone长度
        else:
            print("没有符合条件的序列，跳过综合示意图生成。")

    print("正在生成随机100条综合示意图，包含所有符合条件的序列...")
    read_lengths = get_read_lengths(input_file)
    combined_alignments2 = {read_id: merged_alignments[read_id] for read_id in sample_reads if
                           read_id in merged_alignments}
    plot_all_alignments(combined_alignments2, read_lengths, plot_dir, vector_ref_lengths, backbone_ref_lengths,
                        itr_flipflop_types, "random_100", "Random 100 Sequences (R100S) ")  # 传递Backbone长度

    # 获取匹配上Vector的reads
    # vector_qualifying_reads = get_reads_matched_to_vector(vector_alignments, vector_ref_lengths)

    # 获取所有reads的长度
    # read_lengths = get_read_lengths(input_file)

    # 获取vector和backbone参考序列的总长度
    # fasta_len = sum(vector_ref_lengths.values()) + sum(ITR_ref_lengths.values())*2
    # 绘制匹配上Vector的reads的长度分布图
    # plot_length_distribution(vector_qualifying_reads, read_lengths, plot_dir, "vector_matched_length_distribution", fasta_len)

    # 清理中间文件（可选）
    # os.remove(fastq_fasta)
    # os.remove(itr_blast_out)
    # os.remove(vector_blast_out)
    # os.remove(backbone_blast_out)  # 可选：删除Backbone BLAST输出
    generate_statistics(input_file, merged_alignments, vector_alignments, backbone_blast_out, plot_dir)
    create_reference_sequences("ITR_flipflop.fasta", "references/vector.fasta", "output_references")
    fastq_files = [
        "flip-minus-flip.fastq", "flop-minus-flip.fastq", "flip-minus-flop.fastq", "flop-minus-flop.fastq",
        "flip-plus-flip.fastq", "flop-plus-flip.fastq", "flop-plus-flop.fastq", "flip-plus-flop.fastq"
    ]
    reference_files = [
        "output_references/Left-flip_reverse_Right-flip.fasta",
        "output_references/Left-flop_reverse_Right-flip.fasta",
        "output_references/Left-flip_reverse_Right-flop.fasta", "output_references/Left-flop_reverse_Right-flop.fasta",
        "output_references/Left-flip_forward_Right-flip.fasta", "output_references/Left-flop_forward_Right-flip.fasta",
        "output_references/Left-flop_forward_Right-flop.fasta", "output_references/Left-flip_forward_Right-flop.fasta"

    ]
    output_dir = "bam_output"

    process_fastq_files(fastq_files, reference_files, output_dir)


    print("All done.")


if __name__ == "__main__":
    main()
