# -*- coding: utf-8 -*-
# @Time    : 2023/04/11 15:16
# @Author  : zzk
# @File    :
# Description:

import sys
import time

import yaml
import os
import xlrd
import glob
import docx
import re


def generate_yaml(name_p, RefInfo_name, DemultiInfo_name, seq_):
    '''
    需要copy过来的TempPl.txt,放在当前目录
    需要输入name_P: 当前处理的序列号
    需要输入RefInfo_name: 参考序列yaml文件名
    需要输入DemultiInfo_name: 分流yaml文件名
    需要输入seq_: “seqkit locate TempP1569.fa -p CATGAAGCGCCACCTTTG” 中的这个序列
    '''
    seq_re = seq_.replace("A", "B").replace("G", "D").replace("T", "A").replace("C", "G").replace("B", "T").replace("D",
                                                                                                                    "C")[
             ::-1]
    ITR_L = 'CCTGCAGGCAGCTGCGCGCTCGCTC'
    ITR_R = 'GAGCGAGCGCGCAGCTGCCTGCAGG'
    ltr_l_seq_r = ''
    backbone = ''
    unique_seq = ""
    '''
    1、L在左，R在右
    2、R在左，L在右
    3、找不到
    '''
    # 1 L在左，R在右
    l_pos = seq_.rfind(ITR_L)
    r_pos = seq_.rfind(ITR_R)
    if l_pos != -1 and l_pos < r_pos:
        ltr_l_seq_r = seq_[l_pos:r_pos+len(ITR_R)]
        backbone = seq_[r_pos+len(ITR_R):]+seq_[:l_pos]
        unique_seq = seq_[l_pos+len(ITR_L):r_pos]
    elif r_pos != -1 and l_pos > r_pos:
        ltr_l_seq_r = seq_[l_pos:]+seq_[:r_pos+len(ITR_R)]
        backbone = seq_[r_pos:l_pos+len(ITR_L)]
        unique_seq = seq_[l_pos+len(ITR_L):]+seq_[:r_pos]
    elif l_pos + r_pos == -2:
        l_pos_re = seq_re.rfind(ITR_L)
        r_pos_re = seq_re.rfind(ITR_R)
        if l_pos_re != -1 and l_pos_re < r_pos_re:
            ltr_l_seq_r = seq_re[l_pos_re:r_pos_re + len(ITR_R)]
            backbone = seq_re[r_pos_re + len(ITR_R):] + seq_re[:l_pos_re]
            unique_seq = seq_re[l_pos_re + len(ITR_L):r_pos_re]
        elif r_pos_re != -1 and l_pos_re > r_pos_re:
            ltr_l_seq_r = seq_re[l_pos_re:] + seq_re[:r_pos_re + len(ITR_R)]
            backbone = seq_re[r_pos_re:l_pos_re + len(ITR_L)]
            unique_seq = seq_re[l_pos_re+ len(ITR_L):] + seq_re[:r_pos_re]
        elif l_pos_re + r_pos_re == -2:
            print(file_word,"怎样都找不到ITR的两端")
    elif l_pos * r_pos < 0:
        print(file_word,"ITR一端找得到一端找不到")

    bb = 'CTGTCAGACCAAGTTTACTC'
    pos_bb = backbone.rfind(bb)
    add_R = backbone[:pos_bb]
    add_L = backbone[pos_bb:]
    unique_seq = unique_seq[len(unique_seq)//8:7*(len(unique_seq)//8)]
    dic_sample_ITR = {
        name_p.replace("-Plasmid",""): {"reference_id": name_p.replace("-Plasmid",""), "left_seqs": '', "right_seqs": '', "BCprimer_F": '', "BClen_F": '',
                 "BCprimer_R": '', "BClen_R": '', "unique_sequence": unique_seq}}
    dic_sample_ITR_ref = {name_p.replace("-Plasmid",""): {"sequence": ltr_l_seq_r}}


    dic_sample_ITR_bb = {
        name_p: {"reference_id": name_p, "left_seqs": '', "right_seqs": '', "BCprimer_F": '', "BClen_F": '',
                 "BCprimer_R": '', "BClen_R": '', "unique_sequence": unique_seq}}
    dic_sample_ITR_bb_ref = {name_p: {"sequence": add_L+ltr_l_seq_r+add_R}}

    with open(f'{DemultiInfo_name}.yaml', 'a', encoding='utf-8') as f:
        yaml.dump(dic_sample_ITR, f, default_style=False, encoding='utf-8', allow_unicode=True)
        # yaml.dump(dic_sample_ITR_bb, f, default_style=False, encoding='utf-8', allow_unicode=True)


    with open(f'{RefInfo_name}.yaml', 'a', encoding='utf-8') as f:
        yaml.dump(dic_sample_ITR_ref, f, default_style=False, encoding='utf-8', allow_unicode=True)
        # yaml.dump(dic_sample_ITR_bb_ref, f, default_style=False, encoding='utf-8', allow_unicode=True)


def read_docx(file_doc):
    '''
    word都是以序列号结尾
    读取docx，生成TempPl.txt
    调用这个返回调整好的序列，直接可以放到Temptxt中
    '''
    # 获取文档对象
    flag = False
    seq = ''
    file = docx.Document(file_doc)
    for para in file.paragraphs:
        para_up = para.text.upper().replace(" ", "").replace("\n", "").replace(" ", "")
        count = len(re.findall("A", para_up)) + len(re.findall("T", para_up)) + len(re.findall("G", para_up)) + len(
            re.findall("C", para_up))
        if len(para_up) == 0:
            continue
        if len(para_up) != 0:
            if count / len(para_up) > 0.95:
                seq += para_up
                flag = True
            if count / len(para_up) < 0.6 and flag == True:
                if len(seq) < 1000:
                    flag = False
                    seq = ''
                else:
                    break
    seq_expect = ""
    for i in seq.upper():
        if i in "ATGC":
            seq_expect+=i
    print(len(seq_expect))
    return seq_expect



if __name__ == '__main__':
    all_file_doc = glob.glob("*.docx")  # 取到所有的word文件
    if len(all_file_doc) == 0:
        print("当前目录没有word文件，程序三秒后关闭！")
        time.sleep(3)
        sys.exit()

    print("如果多次运行处理同一批数据，需要将之前生成的yaml删除，原因是数据时依次不覆盖追加！")
    RefInfo_name = 'refinfo'
    # DemultiInfo_name = input("输入DemultiInfo的yaml名：\n")
    DemultiInfo_name = 'Demuinfo'
    # pattern = r'(\S+)\s+(\S+)'
    pattern = r'[P,S]\d{3,4}'
    for file_word in all_file_doc:
        name_x = re.findall(pattern,file_word)
        if len(name_x)!=1:
            print(file_word,"名字中含有两种序号！")
            break
        name_rest = file_word.split("{}".format(name_x[0]))
        print(name_x)
        longest_str = max(name_rest, key=len)
        name_p = (longest_str.strip(".docx")+name_x[0]).replace(" ","").replace("(","-").replace(")","-")
        print(name_p)
        sequence = read_docx(file_word)
        generate_yaml(name_p=name_p, RefInfo_name=RefInfo_name, DemultiInfo_name=DemultiInfo_name, seq_=sequence)
    print("over............")

