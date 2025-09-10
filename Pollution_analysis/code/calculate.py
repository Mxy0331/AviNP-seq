import pandas as pd

# 预定义染色体顺序
CHROM_ORDER = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'
]

# hg38染色体长度（单位：碱基对）
CHROM_LENGTHS = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chrX': 156040895,
    'chrY': 57227415,
}

def calculate_rpmb(input_file, output_file, sheet_name=0):
    # 读取文件
    df = pd.read_excel(input_file, sheet_name=sheet_name)
    
    # 提取主染色体名称
    df['main_chrom'] = df['CHROM'].str.split('_').str[0]
    
    # 统计各染色体行数
    counts = df['main_chrom'].value_counts().to_dict()
    
    # 计算每百万碱基reads数
    results = []
    for chrom in CHROM_ORDER:
        count = counts.get(chrom, 0)
        length = CHROM_LENGTHS.get(chrom, 1)  # 避免除以0错误
        rpmb = count / (length / 1e6)
        results.append(round(rpmb, 2))  # 保留两位小数
        
    # 写入结果文件
    with open(output_file, 'w') as f:
        f.write(f"{CHROM_ORDER}\n")
        f.write(f"{results}\n")

# 使用示例
calculate_rpmb('detailed_report.xlsx', 'output.txt')