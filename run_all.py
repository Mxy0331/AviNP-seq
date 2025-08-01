import os
import subprocess

# def run_demul():
#     print("=== Step 1: Demultiplexing using demul-AAV.py ===")
#     subprocess.run(["python", "demul-AAV.py"], check=True)
#     print("===> Finished demultiplexing.\n")

def run_demul():
    print("=== Step 1: Demultiplexing using demul-AAV.py ===")
    subprocess.run(["python", "demul-AAV.py"], check=True)
    print("===> Finished demultiplexing.\n")

    # 确保 Demultiplexed 目录存在
    os.makedirs("Demultiplexed", exist_ok=True)

    # 将 demul 输出的 fastq 移动或复制到 Demultiplexed
    for file in os.listdir("output"):
        if file.endswith(".fastq"):
            src = os.path.join("output", file)
            dst = os.path.join("Demultiplexed", file)
            # 复制或移动均可，根据需求选择：
            # shutil.copy(src, dst)     # 复制
            os.rename(src, dst)         # 移动


def run_main():
    print("=== Step 2: Running main.py on Demultiplexed fastq files ===")
    os.makedirs("Demultiplexed", exist_ok=True)
    os.makedirs("output", exist_ok=True)
    
    # 将 Demultiplexed 中的 fastq 移动到当前目录，供 main.py 使用
    # 或者根据需要修改 main.py 逻辑来读取 Demultiplexed 中的文件

    subprocess.run(["python", "main.py"], check=True)
    print("===> Finished structure alignment and classification.\n")

def main():
    run_demul()
    run_main()

if __name__ == "__main__":
    main()
