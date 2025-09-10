import pandas as pd
from Bio import SeqIO

# 读取Excel文件中的ID列
def read_ids_from_excel(excel_file):
    # 使用pandas读取excel文件中的ID列
    df = pd.read_excel(excel_file)
    ids = df['ID'].tolist()  # 假设ID列名为'ID'
    return set(ids)  # 使用set去重

# 筛选出nomatch.fastq中与ID匹配的序列
def filter_fastq_by_ids(fastq_file, ids, output_file):
    # 打开输出文件以写入匹配的序列
    with open(output_file, 'w') as out_handle:
        # 使用SeqIO读取fastq文件
        with open(fastq_file, 'r') as in_handle:
            for record in SeqIO.parse(in_handle, 'fastq'):
                # 检查序列ID是否包含在ids列表中的ID
                if any(id_part in record.id for id_part in ids):
                    SeqIO.write(record, out_handle, 'fastq')

def main():
    excel_file = 'detailed_report.xlsx'  # 你的Excel文件路径
    fastq_file = 'nomatch.fastq'  # 你的FastQ文件路径
    output_file = 'align_to_hg.fastq'  # 输出文件的路径

    # 读取ID
    ids = read_ids_from_excel(excel_file)
    
    # 筛选并保存匹配的序列
    filter_fastq_by_ids(fastq_file, ids, output_file)
    print(f"筛选后的序列已保存到{output_file}")

if __name__ == '__main__':
    main()
