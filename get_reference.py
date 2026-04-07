#!/usr/bin/env python3
"""
参考序列提取脚本
用于从FASTA文件中提取指定染色体区域的序列
"""

import sys
import os


def extract_reference(fasta_file, chromosome, start, end):
    """
    从FASTA文件中提取指定染色体区域的序列
    
    Args:
        fasta_file: FASTA文件路径
        chromosome: 染色体名称
        start: 起始位置（1-based）
        end: 结束位置（1-based）
    
    Returns:
        提取的序列字符串
    """
    sequence = ""
    in_target_chromosome = False
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 解析染色体名称
                header = line[1:].split()[0]
                # 处理不同的染色体命名格式
                # 移除可能的'chr'前缀进行比较
                compare_chromosome = chromosome.replace('chr', '')
                header_chromosome = header.split()[0].replace('chr', '')
                if header_chromosome == compare_chromosome:
                    in_target_chromosome = True
                else:
                    in_target_chromosome = False
            elif in_target_chromosome:
                sequence += line
    
    # 提取指定区域（注意：Python是0-based，而位置是1-based）
    if start > 0 and end <= len(sequence):
        return sequence[start-1:end]
    else:
        return ""


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("用法: python get_reference.py <fasta_file> <chromosome> <start> <end>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    chromosome = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    
    if not os.path.exists(fasta_file):
        print(f"错误: 找不到文件 {fasta_file}")
        sys.exit(1)
    
    sequence = extract_reference(fasta_file, chromosome, start, end)
    print(sequence)
