#!/usr/bin/env python3
"""
解析AB1文件名中的染色体位置，与参考基因组比对，绘制指定位置左右25bp的峰图
使用bwa进行比对
"""

import sys
import re
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from Bio import SeqIO
import tempfile
import shutil

DEFAULT_REF_FILE = os.getenv("DEFAULT_REF_FILE", "/data/reference/hg19/hs37d5.fa")
GET_REFERENCE_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "get_reference.py")
BASE_FONT_FAMILY = "DejaVu Sans Mono"
BASE_FONT_SIZE = 15
LABEL_FONT_SIZE = 13


def parse_position_from_filename(filename):
    """
    从文件名中解析位置信息
    只返回位置，不返回染色体（染色体从bwa比对结果获取）
    """
    # 正则表达式匹配位置（按优先级排序）
    patterns = [
        r'(chr[\w]+):(\d+)',           # chr7:140453136 (最高优先级)
        r'(chr[\w]+)-(\d+)',           # chr7-140453136
        r'(chr[\w]+)_(\d+)',           # chr7_140453136
        r'(\d+):(\d+)-(\d+)',         # 7:140453136-140453236
        r'([\w]+):(\d+)-(\d+)',        # 其他格式
        r'-(\d{4,9})_',               # -139409940_ (前面有连字符，后面有下划线)
        r'_(\d{4,9})_',               # _139409940_ (4-9位数字，前后有下划线)
        r'([A-Z0-9]+)_(\d{4,9})_', # 基因名_位置_ (如 NOTCH1_139409940_)
        r'(\w+)-(\d+)-',               # 基因名-位置- (如 12SrRNA-1095-F) (最低优先级)
    ]
    
    for pattern in patterns:
        match = re.search(pattern, filename)
        if match:
            if len(match.groups()) == 2:
                chrom, pos = match.groups()
                return int(pos)
            elif len(match.groups()) == 3:
                chrom, start, end = match.groups()
                # 返回中间位置
                mid_pos = (int(start) + int(end)) // 2
                return mid_pos
    
    return None


def get_reference_sequence(chrom, pos, window_size=50, ref_file=None):
    """
    获取参考基因组序列
    """
    # 尝试从get_reference.py获取参考序列
    if ref_file is None:
        # 使用默认参考文件路径
        ref_file = DEFAULT_REF_FILE
    
    if not os.path.exists(ref_file):
        print(f"错误: 参考文件 {ref_file} 不存在")
        return None
    
    # 计算实际的起始和结束位置
    start = pos - window_size
    end = pos + window_size
    print(f"获取参考序列: {chrom}:{start}-{end}")
    
    try:
        result = subprocess.run(
            [sys.executable, GET_REFERENCE_SCRIPT, ref_file, chrom, str(max(1, start)), str(end)],
            capture_output=True,
            text=True,
            check=True
        )
        sequence = result.stdout.strip()
        print(f"获取到的参考序列: {sequence}")
        print(f"参考序列长度: {len(sequence)}")
        return sequence
    except Exception as e:
        print(f"获取参考序列失败: {e}")
        return None


def align_with_bwa(sample_seq, ref_file, sample_name="sample"):
    """
    使用bwa进行比对
    
    参数:
        sample_seq: 样本序列
        ref_file: 参考基因组文件路径
        sample_name: 样本名称
    
    返回:
        chrom: 染色体名称
        pos: 比对位置
        cigar: CIGAR字符串
        mapq: 比对质量分数
    """
    # 创建临时目录
    tmpdir = tempfile.mkdtemp()
    
    try:
        # 写入样本序列到FASTA文件
        sample_fasta = os.path.join(tmpdir, "sample.fasta")
        with open(sample_fasta, 'w') as f:
            f.write(f">{sample_name}\n{sample_seq}\n")
        
        # 检查bwa索引是否存在
        ref_index = ref_file
        if not os.path.exists(ref_file + ".bwt"):
            print(f"创建bwa索引: {ref_file}")
            subprocess.run(['bwa', 'index', ref_file], check=True, capture_output=True)
        
        # 使用bwa mem进行比对
        print("使用bwa mem进行比对...")
        result = subprocess.run(
            ['bwa', 'mem', '-T', '30', ref_file, sample_fasta],
            capture_output=True,
            text=True,
            check=True
        )
        
        # 解析SAM输出
        sam_lines = result.stdout.strip().split('\n')
        for line in sam_lines:
            if line.startswith('@'):
                continue  # 跳过header行
            
            fields = line.split('\t')
            if len(fields) < 11:
                continue
            
            flag = int(fields[1])
            if flag & 0x4:  # 未比对
                continue
            
            chrom = fields[2]
            pos = int(fields[3]) - 1  # 转换为0-based
            cigar = fields[5]
            mapq = int(fields[4])
            
            print(f"比对结果: {chrom}:{pos}, CIGAR: {cigar}, MAPQ: {mapq}")
            
            return chrom, pos, cigar, mapq
        
        print("警告: 未找到比对结果")
        return None, None, None, None
    
    except Exception as e:
        print(f"bwa比对失败: {e}")
        return None, None, None, None
    
    finally:
        # 清理临时目录
        shutil.rmtree(tmpdir, ignore_errors=True)


def calculate_target_index(bwa_pos, target_pos, cigar, sample_seq, ref_seq):
    """
    根据bwa比对结果计算目标位置在样本序列中的对应索引
    
    参数:
        bwa_pos: bwa比对的起始位置（0-based）
        target_pos: 目标位置（0-based）
        cigar: CIGAR字符串
        sample_seq: 样本序列
        ref_seq: 参考序列
    
    返回:
        target_idx: 目标位置在样本序列中的索引
    """
    try:
        # 解析CIGAR字符串
        import re
        cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        
        # 转换为0-based
        ref_pos = bwa_pos  # bwa已经返回的是0-based
        target_pos_0based = target_pos - 1  # 目标位置是1-based
        sample_pos = 0
        
        print(f"bwa比对起始位置 (0-based): {bwa_pos}")
        print(f"目标位置 (1-based): {target_pos}")
        print(f"目标位置 (0-based): {target_pos_0based}")
        print(f"CIGAR: {cigar}")
        
        # 计算参考序列中的偏移（0-based）
        ref_offset = target_pos_0based - ref_pos
        print(f"参考序列中的偏移 (0-based): {ref_offset}")
        
        # 遍历CIGAR操作，计算样本序列中的对应位置
        for length, op in cigar_ops:
            length = int(length)
            
            if op in ['M', '=', 'X']:  # 匹配/替换
                if ref_pos <= target_pos_0based < ref_pos + length:
                    # 目标位置在这个匹配区域内
                    offset = target_pos_0based - ref_pos
                    # 软剪切的碱基不计入参考序列，但计入样本序列
                    # 所以样本序列索引 = 之前的样本位置 + 偏移
                    target_idx = sample_pos + offset
                    # 验证计算结果
                    print(f"目标位置在匹配区域内，样本序列索引: {target_idx}")
                    print(f"样本序列中对应位置的碱基: {sample_seq[target_idx]}")
                    # 显示目标位置周围的碱基
                    start = max(0, target_idx - 10)
                    end = min(len(sample_seq), target_idx + 11)
                    print(f"目标位置周围的碱基: {sample_seq[start:end]}")
                    print(f"目标位置周围的碱基索引: {start} to {end-1}")
                    return target_idx
                ref_pos += length
                sample_pos += length
            elif op == 'D':  # 缺失（参考序列有，样本序列无）
                if ref_pos <= target_pos_0based < ref_pos + length:
                    # 目标位置在缺失区域内
                    print(f"目标位置在缺失区域内，返回None")
                    return None
                ref_pos += length
            elif op == 'I':  # 插入（样本序列有，参考序列无）
                sample_pos += length
            elif op == 'S':  # 软剪切（样本序列有，参考序列无）
                # 软剪切的碱基不计入参考序列，但计入样本序列
                sample_pos += length
            elif op == 'H':  # 硬剪切（样本序列无，参考序列无）
                pass
            elif op == 'P':  # 填充（样本序列无，参考序列无）
                pass
        
        print("目标位置不在比对区域内，返回None")
        return None
    except Exception as e:
        print(f"计算目标位置失败: {e}")
        return None


def detect_heterozygous(traces, base_locs, base_idx, sample_base):
    """
    检测指定位置是否为杂合
    
    参数:
        traces: 峰图数据字典 {'A': array, 'T': array, 'C': array, 'G': array}
        base_locs: 碱基位置数组
        base_idx: 碱基索引
        sample_base: 样本碱基
    
    返回:
        dict: {
            'is_heterozygous': bool,
            'secondary_base': str,  # 次峰碱基
            'secondary_ratio': float,  # 次峰比例
            'peak_heights': dict  # 各碱基峰高
        }
    """
    try:
        loc = base_locs[base_idx]
        window = 3
        
        peak_heights = {}
        for nt in ['A', 'T', 'C', 'G']:
            trace = traces[nt]
            region = trace[max(0, loc-window):min(len(trace), loc+window+1)]
            peak_heights[nt] = float(np.max(region))
        
        sorted_peaks = sorted(peak_heights.items(), key=lambda x: x[1], reverse=True)
        primary = sorted_peaks[0]
        secondary = sorted_peaks[1]
        
        if primary[1] > 0:
            secondary_ratio = secondary[1] / primary[1]
        else:
            secondary_ratio = 0.0
        
        is_heterozygous = secondary_ratio > 0.25
        
        return {
            'is_heterozygous': is_heterozygous,
            'secondary_base': secondary[0],
            'secondary_ratio': secondary_ratio,
            'peak_heights': peak_heights
        }
    except Exception as e:
        print(f"检测杂合失败: {e}")
        return {
            'is_heterozygous': False,
            'secondary_base': 'N',
            'secondary_ratio': 0.0,
            'peak_heights': {}
        }


def plot_genome_alignment(ab1_file, output_file='genome_alignment.png', 
                        window_size=10, ref_file=None, chrom=None, pos=None,
                        variant_output_file=None):
    """
    与基因组比对并绘制指定位置左右10bp的峰图
    使用bwa进行比对
    
    参数:
        ab1_file: AB1文件路径
        output_file: 输出文件路径
        window_size: 显示窗口大小（默认10bp）
        ref_file: 参考基因组文件路径
        chrom: 手动指定的染色体（可选，如果未指定则使用bwa比对结果）
        pos: 手动指定的位置（可选，如果未指定则使用bwa比对结果）
        variant_output_file: 变异结果输出文件路径（可选，默认为output_file同名txt文件）
    """
    print(f"正在处理文件: {ab1_file}")
    
    # 设置默认参考基因组文件路径
    if ref_file is None:
        ref_file = DEFAULT_REF_FILE
    
    print(f"使用参考基因组: {ref_file}")
    
    # 读取AB1文件
    try:
        record = SeqIO.read(ab1_file, 'abi')
    except Exception as e:
        print(f"读取AB1文件失败: {e}")
        return None
    
    # 获取样本序列
    sample_seq = str(record.seq)
    print(f"样本序列长度: {len(sample_seq)}")
    
    # 保存原始样本序列（用于峰图显示）
    original_sample_seq = sample_seq
    
    # 检查是否为反向测序（R端），如果是则反向互补
    from Bio.Seq import Seq
    # 检查文件名中是否包含-F或-R（正向/反向测序）
    # 通过样本名称中的-F或-R来识别，如 _GA-F_ 或 _GA-R_
    import re
    # 匹配 _XX-F_ 或 _XX-R_ 模式，支持 NF/NR 格式
    # F端: -F_ 或 NF_ 或 -F. 或 NF.
    # R端: -R_ 或 NR_ 或 -R. 或 NR.
    f_match = re.search(r'(NF|-F)([_.-]|$)', ab1_file)
    r_match = re.search(r'(NR|-R)([_.-]|$)', ab1_file)
    is_reverse = bool(r_match) and not bool(f_match)
    if is_reverse:
        print("检测到反向测序（R端），进行反向互补")
        sample_seq = str(Seq(sample_seq).reverse_complement())
        print(f"反向互补后的样本序列长度: {len(sample_seq)}")
    
    # 使用bwa进行比对，获取染色体信息
    sample_name = os.path.basename(ab1_file).replace('.ab1', '')
    bwa_chrom, bwa_pos, cigar, mapq = align_with_bwa(sample_seq, ref_file, sample_name)
    
    if bwa_chrom is None:
        print("bwa比对失败，无法确定染色体")
        return None
    
    # 使用bwa获取的染色体信息
    if chrom is None:
        chrom = bwa_chrom
        print(f"使用bwa获取的染色体: {chrom}")
    else:
        print(f"使用手动指定的染色体: {chrom}")
    
    # 解析CIGAR字符串，确定匹配的长度
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    # 计算参考序列长度（包括缺失操作D）
    ref_length = 0
    for length, op in cigar_ops:
        if op in ['M', '=', 'X', 'D']:  # 匹配或缺失操作
            ref_length += int(length)
    
    # 使用手动输入的位置信息
    if pos is None:
        # 先尝试从文件名中解析位置
        filename = os.path.basename(ab1_file)
        pos_from_filename = parse_position_from_filename(filename)
        
        if pos_from_filename is not None:
            pos = pos_from_filename
            print(f"从文件名解析位置: {pos}")
        else:
            # 如果文件名中没有位置信息，使用BWA比对结果的中心位置
            # bwa_pos是0-based，ref_length是参考序列长度
            # 中心位置 = bwa_pos + ref_length / 2
            # 转换为1-based坐标
            center_pos = bwa_pos + ref_length // 2
            pos = center_pos
            print(f"使用bwa比对结果的中心位置: {pos} (起点: {bwa_pos}, 终点: {bwa_pos + ref_length})")
    else:
        print(f"使用手动指定的位置: {pos}")
    
    print(f"最终使用的染色体位置: {chrom}:{pos}")
    
    # 获取参考序列（使用bwa比对位置和参考序列长度）
    # 注意：bwa返回的是0-based坐标，但get_reference.py使用1-based坐标
    # 所以需要将0-based坐标转换为1-based坐标
    ref_start_0based = bwa_pos
    ref_end_0based = bwa_pos + ref_length
    
    # 转换为1-based坐标用于get_reference.py
    ref_start_1based = max(1, ref_start_0based + 1)
    ref_end_1based = ref_end_0based + 1
    
    # 处理染色体名称转换（chrM -> MT）
    ref_chrom = chrom
    if chrom.lower() == 'chrm':
        ref_chrom = 'MT'
    
    # 获取参考序列
    print(f"获取参考序列: {ref_chrom}:{ref_start_1based}-{ref_end_1based}")
    try:
        result = subprocess.run(
            [sys.executable, GET_REFERENCE_SCRIPT, ref_file, ref_chrom, str(ref_start_1based), str(ref_end_1based)],
            capture_output=True,
            text=True,
            check=True
        )
        ref_seq = result.stdout.strip()
        print(f"获取到的参考序列: {ref_seq}")
        print(f"参考序列长度: {len(ref_seq)}")
    except Exception as e:
        print(f"获取参考序列失败: {e}")
        print("无法获取参考序列，使用样本序列作为参考")
        ref_seq = sample_seq
    
    print(f"参考序列长度: {len(ref_seq)}")
    
    # 直接使用bwa的CIGAR字符串来提取匹配的部分，而不是重新比对
    # 解析CIGAR字符串，确定样本序列和参考序列的匹配区域
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    sample_start = 0
    ref_start = 0
    sample_matched = []
    ref_matched = []
    variants = []
    
    current_sample_pos = 0
    current_ref_pos = 0
    
    for length, op in cigar_ops:
        length = int(length)
        
        if op in ['M', '=', 'X']:  # 匹配/替换
            for i in range(length):
                sample_base = sample_seq[current_sample_pos + i]
                ref_base = ref_seq[current_ref_pos + i]
                sample_matched.append(sample_base)
                ref_matched.append(ref_base)
                
                # 检查变异 (记录原始样本序列位置)
                if sample_base != ref_base:
                    variants.append((current_sample_pos + i, 'mutation', ref_base, sample_base))
            
            current_sample_pos += length
            current_ref_pos += length
        
        elif op == 'D':  # 缺失（参考序列有，样本序列无）
            for i in range(length):
                ref_base = ref_seq[current_ref_pos + i]
                sample_matched.append('-')
                ref_matched.append(ref_base)
                variants.append((current_sample_pos + i - 1 if current_sample_pos > 0 else 0, 'deletion', ref_base, '-'))
            current_ref_pos += length
        
        elif op == 'I':  # 插入（样本序列有，参考序列无）
            for i in range(length):
                sample_base = sample_seq[current_sample_pos + i]
                sample_matched.append(sample_base)
                ref_matched.append('-')
                variants.append((current_sample_pos + i, 'insertion', '-', sample_base))
            current_sample_pos += length
        
        elif op == 'S':  # 软剪切（样本序列有，参考序列无）
            # 软剪切不参与比对，跳过
            current_sample_pos += length
        
        elif op == 'H':  # 硬剪切（样本序列无，参考序列无）
            pass
        
        elif op == 'P':  # 填充（样本序列无，参考序列无）
            pass
    
    aligned_sample = ''.join(sample_matched)
    aligned_ref = ''.join(ref_matched)
    sample_end = current_sample_pos
    ref_end = current_ref_pos
    
    print(f"样本序列比对区域: {sample_start}-{sample_end}")
    print(f"参考序列比对区域: {ref_start}-{ref_end}")
    print(f"比对区域长度: {sample_end - sample_start}")
    print(f"比对后长度: {len(aligned_sample)}")
    print(f"发现 {len(variants)} 个变异")
    
    # 分析比对结果
    if len(aligned_sample) > 0:
        # 计算比对得分（只计算匹配的碱基）
        match_count = sum(1 for a, b in zip(aligned_sample, aligned_ref) if a == b and a != '-')
        alignment_score = match_count / len(aligned_sample) * 100
        print(f"比对相似度: {alignment_score:.2f}%")
    
    # 获取峰图数据
    traces = {
        'A': np.array(record.annotations['abif_raw']['DATA10']),
        'T': np.array(record.annotations['abif_raw']['DATA11']),
        'C': np.array(record.annotations['abif_raw']['DATA12']),
        'G': np.array(record.annotations['abif_raw']['DATA9'])
    }
    
    # 获取碱基位置
    base_locs = np.array(record.annotations['abif_raw']['PLOC1'])
    
    # 对检测到的变异进行杂合检测
    enhanced_variants = []
    for var in variants:
        var_pos, var_type, ref_base, var_base = var[:4]
        # 计算变异在原始样本序列中的索引
        # 需要考虑软剪切
        actual_sample_idx = var_pos
        
        # 检测杂合
        if var_type == 'mutation' and actual_sample_idx < len(base_locs):
            het_result = detect_heterozygous(traces, base_locs, actual_sample_idx, var_base)
            enhanced_variants.append({
                'position': var_pos,
                'type': var_type,
                'ref_base': ref_base,
                'var_base': var_base,
                'is_heterozygous': het_result['is_heterozygous'],
                'secondary_base': het_result['secondary_base'],
                'secondary_ratio': het_result['secondary_ratio'],
                'peak_heights': het_result['peak_heights']
            })
            if het_result['is_heterozygous']:
                print(f"变异位置 {var_pos}: {ref_base}->{var_base} (杂合, 次峰: {het_result['secondary_base']}, 比例: {het_result['secondary_ratio']:.2f})")
            else:
                print(f"变异位置 {var_pos}: {ref_base}->{var_base} (纯合)")
        else:
            enhanced_variants.append({
                'position': var_pos,
                'type': var_type,
                'ref_base': ref_base,
                'var_base': var_base,
                'is_heterozygous': False,
                'secondary_base': 'N',
                'secondary_ratio': 0.0,
                'peak_heights': {}
            })
    
    variants = enhanced_variants
    
    # 输出变异结果到txt文件
    if variant_output_file is None and output_file:
        variant_output_file = os.path.splitext(output_file)[0] + '_variants.txt'
    
    if variants and variant_output_file:
        with open(variant_output_file, 'w', encoding='utf-8') as f:
            f.write(f"# 变异检测结果\n")
            f.write(f"# 文件: {ab1_file}\n")
            f.write(f"# 染色体: {chrom}\n")
            f.write(f"# 比对位置: {bwa_pos}\n")
            f.write(f"# CIGAR: {cigar}\n")
            f.write(f"# 变异数量: {len(variants)}\n")
            f.write(f"#\n")
            f.write(f"样本位置\t类型\t参考碱基\t样本碱基\t纯杂合状态\t次峰碱基\t次峰比例\t峰高(A)\t峰高(T)\t峰高(C)\t峰高(G)\t基因组位置\n")
            
            for var in variants:
                var_pos = var['position']
                var_type = var['type']
                ref_base = var['ref_base']
                var_base = var['var_base']
                is_het = var['is_heterozygous']
                sec_base = var['secondary_base']
                sec_ratio = var['secondary_ratio']
                peak_heights = var['peak_heights']
                
                het_status = "杂合" if is_het else "纯合"
                
                peak_a = peak_heights.get('A', 0)
                peak_t = peak_heights.get('T', 0)
                peak_c = peak_heights.get('C', 0)
                peak_g = peak_heights.get('G', 0)
                
                # 计算基因组位置
                if bwa_pos and cigar:
                    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
                    ref_offset = 0
                    sample_offset = 0
                    for length, op in cigar_ops:
                        length = int(length)
                        if op in ['M', '=', 'X']:
                            if sample_offset + length > var_pos:
                                genome_pos = bwa_pos + ref_offset + (var_pos - sample_offset)
                                break
                            sample_offset += length
                            ref_offset += length
                        elif op == 'D':
                            ref_offset += length
                        elif op == 'I':
                            if sample_offset + length > var_pos:
                                genome_pos = bwa_pos + ref_offset
                                break
                            sample_offset += length
                        elif op == 'S':
                            if sample_offset + length > var_pos:
                                genome_pos = bwa_pos + ref_offset
                                break
                            sample_offset += length
                    else:
                        genome_pos = bwa_pos + ref_offset
                else:
                    genome_pos = "N/A"
                
                f.write(f"{var_pos}\t{var_type}\t{ref_base}\t{var_base}\t{het_status}\t{sec_base}\t{sec_ratio:.2f}\t{peak_a:.0f}\t{peak_t:.0f}\t{peak_c:.0f}\t{peak_g:.0f}\t{genome_pos}\n")
        
        print(f"变异结果已保存到: {variant_output_file}")
    
    # 确定目标位置在样本序列中的索引
    # 使用bwa比对结果计算目标位置
    target_idx = None
    ref_seq_start = pos - len(sample_seq)
    target_offset_in_ref = pos - ref_seq_start
    
    if cigar and bwa_pos is not None:
        # 解析CIGAR字符串，计算目标位置在样本序列中的对应位置
        target_idx = calculate_target_index(bwa_pos, pos, cigar, sample_seq, ref_seq)
    
    if target_idx is None:
        # 备用计算方法
        # 根据参考序列中的目标位置计算样本序列中的对应位置
        # 参考序列的起始位置是 pos - len(sample_seq)
        ref_seq_start = pos - len(sample_seq)
        # 目标位置在参考序列中的偏移量
        target_offset_in_ref = pos - ref_seq_start
        # 根据比对结果，计算目标位置在样本序列中的索引
        target_idx = sample_start + (target_offset_in_ref - ref_start)
    
    # 如果是反向测序，需要将目标位置索引转换回原始序列
    if is_reverse:
        # 原始序列的第i个碱基对应反向互补序列的第(N-1-i)个碱基
        # 所以反向互补序列的第target_idx个碱基对应原始序列的第(N-1-target_idx)个碱基
        target_idx = len(sample_seq) - 1 - target_idx
        print(f"反向测序：转换目标位置索引到原始序列: {target_idx}")
    
    # 确保索引在有效范围内
    target_idx = min(max(0, target_idx), len(base_locs) - 1)
    
    print(f"参考序列起始位置: {ref_seq_start}")
    print(f"目标位置在参考序列中的偏移: {target_offset_in_ref}")
    print(f"目标位置在样本序列中的索引: {target_idx}")
    
    # 显示目标位置周围的碱基（使用原始序列）
    target_surround_start = max(0, target_idx - 10)
    target_surround_end = min(len(original_sample_seq), target_idx + 11)
    target_surrounding_bases = original_sample_seq[target_surround_start:target_surround_end]
    print(f"目标位置周围的碱基: {target_surrounding_bases}")
    print(f"目标位置周围的碱基索引: {target_surround_start} to {target_surround_end}")
    
    # 确定显示范围（左右10bp）
    start_idx = max(0, target_idx - window_size)
    end_idx = min(len(base_locs), target_idx + window_size + 1)
    
    # 使用原始样本序列进行显示（确保与峰图对应）
    display_sequence = original_sample_seq[start_idx:end_idx]
    display_base_locs = base_locs[start_idx:end_idx]
    
    print(f"显示范围: {start_idx}-{end_idx}")
    print(f"显示序列: {display_sequence}")
    
    # 确定x轴范围
    x_start = display_base_locs[0] - 5
    x_end = display_base_locs[-1] + 5
    x_start = max(0, x_start)
    x_end = min(len(traces['A']), x_end)
    
    # 创建图表
    fig, (ax_seq, ax_trace, ax_ref) = plt.subplots(3, 1, figsize=(16, 9), 
                                          gridspec_kw={'height_ratios': [1, 3, 1]}, 
                                          sharex=True)
    
    # 添加标题和CIGAR/MAPQ信息
    plt.suptitle(f"Chromosome: {chrom}:{pos}", fontsize=14, y=0.98)
    plt.figtext(0.5, 0.93, f"CIGAR: {cigar}, MAPQ: {mapq}", ha='center', fontsize=12, color='gray')
    
    # 颜色配置
    BASE_COLORS = {
        'A': '#00CC00',  # 绿色
        'T': '#FF0000',  # 红色
        'C': '#0000FF',  # 蓝色
        'G': '#000000',  # 黑色
        'N': '#808080',  # 灰色
        '-': '#CCCCCC'   # 灰色（gap）
    }
    
    # 上方：序列显示
    ax_seq.set_ylim(0, 2)
    ax_seq.set_yticks([])
    
    # 隐藏x轴刻度（只保留峰图的刻度）
    ax_seq.set_xticks([])
    ax_seq.set_xticklabels([])
    ax_seq.tick_params(bottom=False, top=False, labelbottom=False, labeltop=False, which='both')
    
    # 隐藏边框
    ax_seq.spines['top'].set_visible(False)
    ax_seq.spines['right'].set_visible(False)
    ax_seq.spines['bottom'].set_visible(False)
    ax_seq.spines['left'].set_visible(False)
    
    # 绘制样本序列（使用彩色文本方式）
    for i, (base, base_loc) in enumerate(zip(display_sequence, display_base_locs)):
        actual_idx = start_idx + i
        color = BASE_COLORS.get(base, '#808080')
        
        # 绘制彩色碱基
        ax_seq.text(base_loc, 1.5, base, ha='center', va='center',
                  color=color, fontsize=BASE_FONT_SIZE, fontweight='bold',
                  fontfamily=BASE_FONT_FAMILY)
    
    # 添加文件名称标签
    file_name = os.path.basename(ab1_file)
    # 提取文件名中的关键部分（如QW25010306）
    # 格式：QW25010306_(4121F)12SrRNA-1095-F_WG250828057135.ab1 -> QW25010306
    import re
    match = re.match(r'^([A-Za-z0-9]+)[_\(]', file_name)
    if match:
        sample_name = match.group(1)
    else:
        # 如果没有匹配到，使用文件名去掉扩展名
        sample_name = os.path.splitext(file_name)[0]
    ax_seq.text(display_base_locs[0] - 20, 1.5, sample_name, ha='right', va='center',
              fontsize=LABEL_FONT_SIZE, fontweight='bold', color='#333333',
              fontfamily=BASE_FONT_FAMILY)
    
    # 下方：参考序列显示
    ax_ref.set_ylim(0, 2)
    ax_ref.set_yticks([])
    
    # 隐藏参考序列轴的刻度
    ax_ref.set_xticks([])
    ax_ref.set_xticklabels([])
    ax_ref.tick_params(bottom=False, labelbottom=False, which='both')
    
    # 隐藏边框
    ax_ref.spines['top'].set_visible(False)
    ax_ref.spines['right'].set_visible(False)
    ax_ref.spines['bottom'].set_visible(False)
    ax_ref.spines['left'].set_visible(False)
    
    # 获取参考序列的对应区域
    # 计算参考序列中对应于样本序列显示区域的位置
    # 样本序列显示区域：[start_idx, end_idx)
    # 样本序列比对区域：[sample_start, sample_end)
    # 参考序列比对区域：[ref_start, ref_end)
    
    # 计算样本序列显示区域相对于比对区域的偏移
    sample_display_offset = start_idx - sample_start
    sample_display_length = end_idx - start_idx
    
    # 获取参考序列的对应区域
    # 根据CIGAR字符串中的软剪切（S）数量，在参考序列左右补充对应数量的碱基
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    # 计算左端和右端的软剪切数量
    left_soft_clipping = 0
    right_soft_clipping = 0
    found_match = False
    
    for length, op in cigar_ops:
        length = int(length)
        if op == 'S' and not found_match:
            # 第一个S操作是左端软剪切
            left_soft_clipping = length
        elif op == 'S':
            # 后续的S操作是右端软剪切
            right_soft_clipping = length
        elif op in ['M', '=', 'X', 'I', 'D']:
            # 遇到匹配或插入/缺失操作，说明已经过了左端软剪切
            found_match = True
    
    print(f"CIGAR: {cigar}")
    print(f"左端软剪切: {left_soft_clipping}, 右端软剪切: {right_soft_clipping}")
    
    # 如果是反向测序，需要将原始序列索引转换为反向互补序列索引
    if is_reverse:
        # 原始序列的第i个碱基对应反向互补序列的第(N-1-i)个碱基
        # 所以原始序列显示区域[start_idx, end_idx)对应反向互补序列显示区域[N-end_idx, N-start_idx)
        seq_len = len(sample_seq)
        rev_start_idx = seq_len - end_idx
        rev_end_idx = seq_len - start_idx
        
        # 计算反向互补序列显示区域相对于比对区域的偏移
        # 注意：比对区域是从left_soft_clipping开始的
        rev_sample_display_offset = rev_start_idx - left_soft_clipping
        rev_sample_display_length = rev_end_idx - rev_start_idx
        
        # 计算参考序列中对应的显示区域（R端）
        # 参考序列区域 = 参考序列起始位置 + 反向互补序列显示偏移
        ref_display_start = ref_start + rev_sample_display_offset
        ref_display_end = ref_start + rev_sample_display_offset + rev_sample_display_length
        
        # 反向互补参考序列以与原始样本序列对应
        ref_display_sequence = ref_seq[ref_display_start:ref_display_end]
        ref_display_sequence = str(Seq(ref_display_sequence).reverse_complement())
        print(f"反向测序：参考序列显示区域 [{ref_display_start}, {ref_display_end})，反向互补")
    else:
        # 正向测序，直接对应，左右补充软剪切数量的碱基
        ref_display_start = ref_start + sample_display_offset - left_soft_clipping
        ref_display_end = ref_start + sample_display_offset + sample_display_length + right_soft_clipping
        
        # 获取参考序列的显示区域
        ref_display_sequence = ref_seq[ref_display_start:ref_display_end]
        print(f"正向测序：参考序列显示区域 [{ref_display_start}, {ref_display_end})，左右补充软剪切碱基")
    
    # 确保参考序列索引在有效范围内
    ref_display_start = max(0, ref_display_start)
    ref_display_end = min(len(ref_seq), ref_display_end)
    
    # 确保参考序列长度与显示区域匹配
    if len(ref_display_sequence) < len(display_sequence):
        # 如果参考序列长度不足，用'-'填充
        ref_display_sequence = ref_display_sequence.ljust(len(display_sequence), '-')
    elif len(ref_display_sequence) > len(display_sequence):
        # 如果参考序列长度过长，截断
        ref_display_sequence = ref_display_sequence[:len(display_sequence)]
    
    # 调试信息
    print(f"样本序列显示区域偏移: {sample_display_offset}")
    print(f"参考序列显示区域起始: {ref_display_start}")
    print(f"参考序列显示区域结束: {ref_display_end}")
    print(f"参考序列显示序列: {ref_display_sequence}")
    
    # 绘制参考序列（使用彩色文本方式）
    for i, (base, base_loc) in enumerate(zip(ref_display_sequence, display_base_locs)):
        color = BASE_COLORS.get(base, '#808080')
        
        # 绘制彩色碱基
        ax_ref.text(base_loc, 1.5, base, ha='center', va='center',
                  color=color, fontsize=BASE_FONT_SIZE, fontweight='bold',
                  fontfamily=BASE_FONT_FAMILY)
    
    # 添加参考序列标签
    ax_ref.text(display_base_locs[0] - 20, 1.5, 'Reference:', ha='right', va='center',
              fontsize=LABEL_FONT_SIZE, fontweight='bold', color='#333333',
              fontfamily=BASE_FONT_FAMILY)
    
    # 下方：峰图显示
    x = np.arange(x_start, x_end)
    
    # 绘制峰图
    for base in ['A', 'T', 'C', 'G']:
        y = traces[base][x_start:x_end]
        ax_trace.plot(x, y, color=BASE_COLORS[base], linewidth=1.5, alpha=0.8, label=base)
    
    # 在中心位置添加红色虚线
    target_in_display = target_idx - start_idx
    if 0 <= target_in_display < len(display_base_locs):
        target_base_loc = display_base_locs[target_in_display]
        ax_seq.axvline(x=target_base_loc, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        ax_trace.axvline(x=target_base_loc, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        ax_ref.axvline(x=target_base_loc, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # 隐藏峰图边框
    ax_trace.spines['top'].set_visible(False)
    ax_trace.spines['right'].set_visible(False)
    ax_trace.spines['bottom'].set_visible(False)
    ax_trace.spines['left'].set_visible(False)
    
    # 底部位置刻度
    tick_positions = []
    tick_labels = []
    for i, base_loc in enumerate(display_base_locs):
        actual_pos = start_idx + i
        if actual_pos % 5 == 0:
            tick_positions.append(base_loc)
            tick_labels.append(str(actual_pos + 1))
    
    # 设置横坐标刻度（显示在峰图顶部）
    ax_trace.set_xticks(tick_positions)
    ax_trace.set_xticklabels(tick_labels, fontsize=10)
    
    # 将刻度和标签移到顶部，强制隐藏底部刻度
    ax_trace.xaxis.set_ticks_position('top')
    ax_trace.xaxis.set_label_position('top')
    ax_trace.tick_params(bottom=False, labelbottom=False, which='both')
    
    # 设置标题和标签
    ax_seq.set_title("", fontsize=14, fontweight='bold', pad=20)
    
    # 取消x轴标签
    # ax_trace.set_xlabel('Sample Position', fontsize=12, fontweight='bold')
    ax_trace.legend(loc='upper right', fontsize=10, framealpha=0.9)
    
    # 调整布局
    plt.subplots_adjust(hspace=0.1, top=0.85, bottom=0.15)
    
    # 保存图片（如果指定了输出文件）
    if output_file:
        plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor='white')
        print(f"基因组比对峰图已保存到: {output_file}")
    plt.close()
    
    # 返回详细的结果数据
    return {
        'chrom': chrom,
        'pos': pos,
        'window_size': window_size,
        'sample_sequence': sample_seq,
        'reference_sequence': ref_seq,
        'variants': variants,
        'alignment_score': alignment_score if 'alignment_score' in locals() else None,
        'cigar': cigar,
        'mapq': mapq,
        'sample_length': len(sample_seq),
        'reference_length': len(ref_seq),
        'is_reverse': is_reverse,
        'bwa_pos': bwa_pos,
        'target_idx': target_idx
    }


def batch_process_ab1_files(ab1_files, output_dir=None, window_size=25, ref_file=None, combined_output=None):
    """
    批量处理多个AB1文件
    
    参数:
        ab1_files: AB1文件列表
        output_dir: 输出目录
        window_size: 显示窗口大小
        ref_file: 参考基因组文件
        combined_output: 合并输出文件路径（可选，默认为output_dir/variants_combined.txt）
    """
    if output_dir is None:
        output_dir = os.getcwd()
    
    os.makedirs(output_dir, exist_ok=True)
    
    if combined_output is None:
        combined_output = os.path.join(output_dir, 'variants_combined.txt')
    
    all_variants = []
    results = []
    for ab1_file in ab1_files:
        if not ab1_file.endswith('.ab1'):
            print(f"跳过非AB1文件: {ab1_file}")
            continue
        
        filename = os.path.basename(ab1_file)
        output_file = os.path.join(output_dir, f"genome_alignment_{os.path.splitext(filename)[0]}.png")
        
        print(f"\n处理文件: {filename}")
        result = plot_genome_alignment(ab1_file, output_file, window_size, ref_file)
        if result:
            results.append(result)
            # 收集变异信息
            for var in result.get('variants', []):
                all_variants.append({
                    'file': ab1_file,
                    'filename': filename,
                    'chrom': result.get('chrom'),
                    'bwa_pos': result.get('bwa_pos'),
                    'cigar': result.get('cigar'),
                    **var
                })
    
    # 输出合并的变异结果
    if all_variants:
        with open(combined_output, 'w', encoding='utf-8') as f:
            f.write(f"# 变异检测结果汇总\n")
            f.write(f"# 处理文件数: {len(results)}\n")
            f.write(f"# 变异总数: {len(all_variants)}\n")
            f.write(f"# 生成时间: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"#\n")
            f.write(f"文件名\t染色体\t样本位置\t类型\t参考碱基\t样本碱基\t纯杂合状态\t次峰碱基\t次峰比例\t峰高(A)\t峰高(T)\t峰高(C)\t峰高(G)\t基因组位置\n")
            
            for var in all_variants:
                filename = var.get('filename', '')
                chrom = var.get('chrom', '')
                var_pos = var.get('position', '')
                var_type = var.get('type', '')
                ref_base = var.get('ref_base', '')
                var_base = var.get('var_base', '')
                is_het = var.get('is_heterozygous', False)
                sec_base = var.get('secondary_base', 'N')
                sec_ratio = var.get('secondary_ratio', 0.0)
                peak_heights = var.get('peak_heights', {})
                bwa_pos = var.get('bwa_pos')
                cigar = var.get('cigar', '')
                
                het_status = "杂合" if is_het else "纯合"
                
                peak_a = peak_heights.get('A', 0)
                peak_t = peak_heights.get('T', 0)
                peak_c = peak_heights.get('C', 0)
                peak_g = peak_heights.get('G', 0)
                
                # 计算基因组位置
                if bwa_pos and cigar:
                    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
                    ref_offset = 0
                    sample_offset = 0
                    genome_pos = bwa_pos
                    for length, op in cigar_ops:
                        length = int(length)
                        if op in ['M', '=', 'X']:
                            if sample_offset + length > var_pos:
                                genome_pos = bwa_pos + ref_offset + (var_pos - sample_offset)
                                break
                            sample_offset += length
                            ref_offset += length
                        elif op == 'D':
                            ref_offset += length
                        elif op == 'I':
                            if sample_offset + length > var_pos:
                                genome_pos = bwa_pos + ref_offset
                                break
                            sample_offset += length
                        elif op == 'S':
                            if sample_offset + length > var_pos:
                                genome_pos = bwa_pos + ref_offset
                                break
                            sample_offset += length
                    else:
                        genome_pos = bwa_pos + ref_offset
                else:
                    genome_pos = "N/A"
                
                f.write(f"{filename}\t{chrom}\t{var_pos}\t{var_type}\t{ref_base}\t{var_base}\t{het_status}\t{sec_base}\t{sec_ratio:.2f}\t{peak_a:.0f}\t{peak_t:.0f}\t{peak_c:.0f}\t{peak_g:.0f}\t{genome_pos}\n")
        
        print(f"\n合并变异结果已保存到: {combined_output}")
    
    print(f"\n共处理 {len(results)} 个AB1文件")
    return results


def main():
    import subprocess
    
    if len(sys.argv) < 2:
        print("用法: python plot_genome_alignment.py <ab1_file> [chrom:pos] [output_file] [window_size] [ref_file]")
        print("")
        print("参数:")
        print("  ab1_file: AB1文件路径")
        print("  chrom:pos: 手动指定染色体位置（可选）")
        print("  output_file: 输出文件路径（可选）")
        print("  window_size: 显示窗口大小（可选，默认10bp）")
        print("  ref_file: 参考基因组文件路径（可选）")
        print("")
        print("示例:")
        print("  python plot_genome_alignment.py sample.ab1")
        print("  python plot_genome_alignment.py sample.ab1 chr6:32006387")
        print("  python plot_genome_alignment.py sample.ab1 chr6:32006387 alignment.png 10")
        print("")
        print("批量处理多个文件:")
        print("  python plot_genome_alignment.py *.ab1")
        sys.exit(1)
    
    # 检查是否为多个文件
    if '*' in sys.argv[1] or os.path.isdir(sys.argv[1]):
        import glob
        if os.path.isdir(sys.argv[1]):
            ab1_files = glob.glob(os.path.join(sys.argv[1], '*.ab1'))
        else:
            ab1_files = glob.glob(sys.argv[1])
        
        # 解析命名参数
        output_dir = 'genome_alignments'
        window_size = 10
        ref_file = None
        combined_output = None
        
        i = 2
        while i < len(sys.argv):
            if sys.argv[i] == '--window_size' and i + 1 < len(sys.argv):
                window_size = int(sys.argv[i + 1])
                i += 2
            elif sys.argv[i] == '--output_dir' and i + 1 < len(sys.argv):
                output_dir = sys.argv[i + 1]
                i += 2
            elif sys.argv[i] == '--ref_file' and i + 1 < len(sys.argv):
                ref_file = sys.argv[i + 1]
                i += 2
            elif sys.argv[i] == '--combined_output' and i + 1 < len(sys.argv):
                combined_output = sys.argv[i + 1]
                i += 2
            elif not sys.argv[i].startswith('--'):
                # 位置参数
                if output_dir == 'genome_alignments':
                    output_dir = sys.argv[i]
                elif window_size == 10:
                    try:
                        window_size = int(sys.argv[i])
                    except ValueError:
                        pass
                i += 1
            else:
                i += 1
        
        batch_process_ab1_files(ab1_files, output_dir, window_size, ref_file, combined_output)
    else:
        ab1_file = sys.argv[1]
        
        # 检查是否手动指定了染色体位置
        chrom = None
        pos = None
        arg_offset = 2
        
        if len(sys.argv) > 2 and ':' in sys.argv[2] and not sys.argv[2].startswith('--'):
            # 手动指定了染色体位置
            chrom_pos = sys.argv[2]
            parts = chrom_pos.split(':')
            if len(parts) == 2:
                chrom = parts[0]
                try:
                    pos = int(parts[1])
                    arg_offset = 3
                except ValueError:
                    print(f"错误: 无效的位置格式: {parts[1]}")
                    sys.exit(1)
        
        # 解析命名参数
        output_file = f'genome_alignment_{os.path.splitext(os.path.basename(ab1_file))[0]}.png'
        window_size = 10
        ref_file = None
        
        i = arg_offset
        while i < len(sys.argv):
            if sys.argv[i] == '--window_size' and i + 1 < len(sys.argv):
                window_size = int(sys.argv[i + 1])
                i += 2
            elif sys.argv[i] == '--output' and i + 1 < len(sys.argv):
                output_file = sys.argv[i + 1]
                i += 2
            elif sys.argv[i] == '--ref_file' and i + 1 < len(sys.argv):
                ref_file = sys.argv[i + 1]
                i += 2
            elif not sys.argv[i].startswith('--'):
                # 位置参数
                if arg_offset == 2 or (arg_offset == 3 and output_file.startswith('genome_alignment_')):
                    output_file = sys.argv[i]
                elif window_size == 10:
                    try:
                        window_size = int(sys.argv[i])
                    except ValueError:
                        pass
                i += 1
            else:
                i += 1
        
        # 调用函数时传递手动指定的染色体位置
        plot_genome_alignment(ab1_file, output_file, window_size, ref_file, chrom, pos)


if __name__ == "__main__":
    main()
