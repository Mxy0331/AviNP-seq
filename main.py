import os
import subprocess
import shutil
import yaml
import docx
import re
import argparse

def extract_sequence_from_docx(docx_file, start_seq, end_seq):
    """
    从docx文件中提取序列并处理，去除空格和回车，转换为大写，截取指定序列
    :param docx_file: .docx文件路径
    :param start_seq: 截取的起始序列
    :param end_seq: 截取的结束序列
    :return: 处理后的序列
    """
    # 读取docx文件
    doc = docx.Document(docx_file)

    # 将所有文本提取出来，并去除空格和回车
    full_text = ""
    for para in doc.paragraphs:
        full_text += para.text.strip().replace("\n", "").replace(" ", "").upper()

    # 查找并提取指定序列之间的部分
    match = re.search(r"{}(.*?){}".format(re.escape(start_seq), re.escape(end_seq)), full_text)

    if match:
        sequence = match.group(1)
        return sequence
    else:
        print("未找到匹配的序列")
        return None

def run_command(command, cwd=None):
    """帮助执行shell命令"""
    print(f"Running command: {command}")
    try:
        subprocess.run(command, check=True, shell=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        exit(1)

def step_1():
    """执行第一个命令：python ITRyaml.py"""
    run_command("python ITRyaml.py")

def step_2():
    """执行第二个命令：python greporeseq1_6/greporeseq.py all -n N112chop.fastq.gz -d Demuinfo.yaml -r refinfo.yaml"""
    # run_command(f"python greporeseq1_6/greporeseq.py all -n {fq}.fastq.gz -d Demuinfo.yaml -r refinfo.yaml")
    run_command(f"python ref_visualize.py")

    # 2.1 执行 igv.py
    run_command("python igv.py")

    # 2.2 执行 plot_2-new.py
    run_command("python plot_2-new.py")

def step_3(word_name):
    """根据word名自动生成文件夹并复制相关文件"""
    # 创建文件夹
    folder_name = f"{word_name}"
    os.makedirs(folder_name, exist_ok=True)

    # 复制需要的文件
    vector_itr_script = "vector_ITR_Blast_Visualize_6_3.py"
    itr_flipflop_fasta = "ITR_flipflop.fasta"
    ref_yaml = "ref.yaml"
    fastq_folder = "Demultiplexed"  # 当前文件夹

    # 将文件复制到新文件夹
    shutil.copy(vector_itr_script, folder_name)
    shutil.copy(itr_flipflop_fasta, folder_name)
    shutil.copy(ref_yaml, folder_name)

    # 复制相应的fastq文件
    fastq_file = os.path.join(fastq_folder, f"{word_name}.fastq")
    if os.path.exists(fastq_file):
        shutil.copy(fastq_file, folder_name)
    else:
        print(f"Warning: {fastq_file} does not exist.")

    # 修改 ref.yaml 中的 vector 路径
    with open(os.path.join(folder_name, ref_yaml), 'r') as file:
        ref_data = yaml.safe_load(file)

    # 从word里面先把所有的序列的空格和回车去掉。所有的字母大写，截取第一个CCTGCAGGCAGCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT
    # 到第一个AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGCTGCCTGCAGG
    # 之间的序列
    start_seq1 = "CCTGCAGGCAGCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT"
    start_seq2 = "CCTGCAGGCAGCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT"
    end_seq1 = "AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGCTGCCTGCAGG"
    end_seq2 = "AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGCTGCCTGCAGG"

    extracted_sequence = extract_sequence_from_docx(f'{word_name}.docx', start_seq1, end_seq1)
    if not extracted_sequence:
        extracted_sequence = extract_sequence_from_docx(f'{word_name}.docx', start_seq2, end_seq2)
        if not extracted_sequence:
            extracted_sequence = extract_sequence_from_docx(f'{word_name}.docx', start_seq1, end_seq2)
            if not extracted_sequence:
                extracted_sequence = extract_sequence_from_docx(f'{word_name}.docx', start_seq2, end_seq1)
                if not extracted_sequence:
                    print(f"Warning: No valid sequence found in {word_name}.docx")
                    return
    # 假设需要修改的字段是 "vector"，你可以根据实际情况调整
    ref_data['vector'] = {'sequence': extracted_sequence}

    with open(os.path.join(folder_name, ref_yaml), 'w') as file:
        yaml.dump(ref_data, file)

    print(f"Files copied and {ref_yaml} updated in folder: {folder_name}")


def step_4(word_name):
    """进入相应文件夹并执行 vector_ITR_Blast_Visualize_6_3.py"""
    folder_name = f"{word_name}"
    command = f"python vector_ITR_Blast_Visualize_6_3.py -fq {word_name}.fastq -ref ref.yaml --itr_flipflop ITR_flipflop.fasta"

    # 进入对应文件夹并执行命令
    run_command(command, cwd=folder_name)

def step_5():
    """执行 merged_excel.py"""
    run_command("python merged_excel.py")

def step_6():
    run_command('xvfb-run -s "-screen 0 1280x1024x24" igv.sh -b igv.txt')
    print("运行完成！请查看结果文件夹。")

def get_word_names_from_current_folder():
    """自动从当前文件夹获取所有的 word_name (假设是 fastq 文件名不带扩展名)"""
    word_names = []
    for file in os.listdir("."):
        if file.endswith(".docx"):
            word_names.append(file.replace(".docx", ""))
    return word_names

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="自动化执行流程")
    # parser.add_argument('-fq', '--fq', required=True, help="输入的 fastq.gz (例如: N112chop)")
    return parser.parse_args()

def main():
    args = parse_arguments()

    # fq = args.fq
    # 步骤 1
    step_1()

    # 步骤 2
    step_2()

    # 获取 word_names 自动提取来自当前文件夹中的 fastq 文件名
    word_names = get_word_names_from_current_folder()

    for word_name in word_names:
        # 步骤 3
        step_3(word_name)

        # 步骤 4
        step_4(word_name)

    # 步骤 5
    step_5()

    step_6()

if __name__ == "__main__":
    main()
