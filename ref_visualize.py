import os
import yaml
import os
import glob
import subprocess
import os
import glob
import subprocess


def batch_generate_bam(
    ref_dir='reference',
    fq_dir='Demultiplexed',
    out_dir='Visualization',
    threads=4
):
    """
    批量将 Demultiplexed 文件夹下的 FASTQ
    比对到 reference 文件夹下的所有 FASTA，
    生成排序后的 BAM 并建立索引。
    """
    os.makedirs(out_dir, exist_ok=True)

    # 合并所有参考 FASTA 到一个文件
    ref_fasta = os.path.join(ref_dir, 'all_refs.fasta')
    with open(ref_fasta, 'w', encoding='utf-8') as outf:
        for fasta in glob.glob(os.path.join(ref_dir, '*.fasta')):
            with open(fasta, 'r', encoding='utf-8') as inf:
                outf.write(inf.read())

    # 索引参考
    subprocess.run(['bwa', 'index', ref_fasta], check=True)
    subprocess.run(['samtools', 'faidx', ref_fasta], check=True)

    # 对每个 FASTQ 文件进行比对、排序和索引
    for fq in glob.glob(os.path.join(fq_dir, '*.fastq')):
        sample = os.path.basename(fq).rsplit('.fastq', 1)[0]
        unsorted_bam = os.path.join(out_dir, f'{sample}.unsorted.bam')
        sorted_bam = os.path.join(out_dir, f'{sample}.sorted.bam')

        # BWA MEM 比对，输出 BAM
        bwa_proc = subprocess.Popen(
            ['bwa', 'mem', '-t', str(threads), ref_fasta, fq],
            stdout=subprocess.PIPE
        )
        with open(unsorted_bam, 'wb') as bam_out:
            subprocess.run(
                ['samtools', 'view', '-b', '-F', '4', '-'],
                stdin=bwa_proc.stdout,
                stdout=bam_out,
                check=True
            )
        bwa_proc.wait()

        # 排序
        subprocess.run(
            ['samtools', 'sort', '-@', str(threads), '-o', sorted_bam, unsorted_bam],
            check=True
        )
        # 建立索引
        subprocess.run(['samtools', 'index', sorted_bam], check=True)

        print(f'Processed {sample}: {sorted_bam}')


def yaml_to_individual_fastas(yaml_path, output_dir, line_width=60):
    """
    读取一个 refinfo YAML 文件，并为其中每条序列生成一个独立的 FASTA 文件。

    :param yaml_path: 输入 YAML 文件路径
    :param output_dir: 输出目录（会自动创建）
    :param line_width: FASTA 序列每行的长度
    """
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 加载 YAML 数据
    with open(yaml_path, 'r', encoding='utf-8') as yf:
        data = yaml.safe_load(yf)

    # 为每个条目写入独立的 FASTA 文件
    for seq_id, entry in data.items():
        seq = entry.get('sequence', '')
        fasta_path = os.path.join(output_dir, f"{seq_id}.fasta")
        with open(fasta_path, 'w', encoding='utf-8') as ff:
            ff.write(f">{seq_id}\n")
            for i in range(0, len(seq), line_width):
                ff.write(seq[i:i+line_width] + "\n")

if __name__ == "__main__":
    yaml_to_individual_fastas(
        yaml_path="refinfo.yaml",
        output_dir="Reference",
        line_width=60
    )
    print("已在 `reference/` 目录下生成独立的 FASTA 文件。")
    batch_generate_bam()
    print("已在 `Visualization/` 目录下生成排序后的 BAM 文件及索引。")

