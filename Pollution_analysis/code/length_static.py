from collections import defaultdict
import pandas as pd

# 读取FASTQ文件并统计序列长度
def count_read_lengths(fastq_file):
    lengths = defaultdict(int)  # 用来统计不同长度的reads数
    with open(fastq_file, 'r') as f:
        while True:
            # 读取4行，每4行代表一个读段（序列）
            header = f.readline()
            if not header:
                break  # 到文件末尾则退出
            sequence = f.readline().strip()  # 第二行是序列
            f.readline()  # 第三行是质量评分，忽略
            f.readline()  # 第四行是质量评分注释，忽略

            # 统计该序列的长度
            seq_len = len(sequence)
            lengths[seq_len] += 1
    
    return lengths

# 根据长度范围将序列进行分类
def categorize_lengths(lengths):
    categorized = defaultdict(int)
    for seq_len, count in lengths.items():
        # 根据要求的长度区间分组
        if seq_len <= 150:
            category = "0-150"
        elif seq_len <= 300:
            category = "151-300"
        elif seq_len <= 450:
            category = "301-450"
        elif seq_len <= 600:
            category = "451-600"
        elif seq_len <= 750:
            category = "601-750"
        elif seq_len <= 900:
            category = "751-900"
        elif seq_len <= 1050:
            category = "901-1050"
        elif seq_len <= 1200:
            category = "1051-1200"
        elif seq_len <= 1350:
            category = "1201-1350"
        elif seq_len <= 1500:
            category = "1351-1500"
        elif seq_len <= 1650:
            category = "1501-1650"
        elif seq_len <= 1800:
            category = "1651-1800"
        elif seq_len <= 1950:
            category = "1801-1950"
        elif seq_len <= 2100:
            category = "1951-2100"
        elif seq_len <= 2250:
            category = "2101-2250"
        elif seq_len <= 2400:
            category = "2251-2400"
        elif seq_len <= 2500:
            category = "2401-2500"
        else:
            category = ">2500"
        
        categorized[category] += count
    
    return categorized

# 生成txt格式的输出
def generate_txt_output(categorized_lengths):
    # 定义长度范围顺序
    length_ranges = [
        "0-150", "151-300", "301-450", "451-600", "601-750", "751-900",
        "901-1050", "1051-1200", "1201-1350", "1351-1500", "1501-1650",
        "1651-1800", "1801-1950", "1951-2100", "2101-2250", "2251-2400", 
        "2401-2500", ">2500"
    ]
    
    # 获取每个区间的对应读取数值
    counts = [categorized_lengths.get(category, 0) for category in length_ranges]
    
    # 保存到txt文件
    with open("reads_length_distribution.txt", "w") as f:
        f.write(f"[{', '.join(f'\'{range}\'' for range in length_ranges)}]\n")
        f.write(f"[{', '.join(str(count) for count in counts)}]\n")

# 主函数
def main(fastq_file):
    lengths = count_read_lengths(fastq_file)
    categorized_lengths = categorize_lengths(lengths)
    generate_txt_output(categorized_lengths)

if __name__ == "__main__":
    fastq_file = 'align_to_hg.fastq'  # 请替换为你的FASTQ文件路径
    main(fastq_file)
