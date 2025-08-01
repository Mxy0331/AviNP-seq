#!/usr/bin/env python3
import os
import pandas as pd


def collect_and_merge_excel_files(source_dir, output_file):
    """
    遍历 source_dir 及其所有子目录，查找所有文件名以 "statistics.xlsx" 结尾的 Excel 文件，
    将它们的内容合并到一个新的 Excel 文件中（output_file）。
    """
    all_data = []  # 用于存储所有文件的数据

    # 遍历当前目录及所有子目录
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            if file.endswith("statistics.xlsx"):
                file_path = os.path.join(root, file)
                try:
                    # 读取每个 Excel 文件的内容
                    df = pd.read_excel(file_path)
                    all_data.append(df)
                    print(f"已读取文件：{file_path}")
                except Exception as e:
                    print(f"读取文件 {file_path} 时出错: {e}")

    if all_data:
        # 将所有数据合并成一个 DataFrame
        merged_df = pd.concat(all_data, ignore_index=True)

        # 将合并后的数据写入新的 Excel 文件
        merged_df.to_excel(output_file, index=False)
        print(f"合并后的文件已保存为：{output_file}")
    else:
        print("没有找到符合条件的文件。")


if __name__ == "__main__":
    # 获取当前工作目录
    current_dir = os.getcwd()
    # 设置输出的合并文件名
    output_excel_file = os.path.join(current_dir, "merged_statistics.xlsx")

    # 合并所有的 statistics.xlsx 文件
    collect_and_merge_excel_files(current_dir, output_excel_file)
