import sys
from pathlib import Path
import numpy as np

_PY_DIR = Path(__file__).resolve().parent
_BUILD_LIB = _PY_DIR / 'build' / f'lib.win-amd64-cpython-{sys.version_info.major}{sys.version_info.minor}'
if _BUILD_LIB.exists():
    sys.path.insert(0, str(_BUILD_LIB))

import fftqspa
from Inner_Code_DNA_Channel_Simulation import DNAChannel
import os
import csv

CODE_PARM = 2  # i for 2^i-ary LDPC, i=1, 2, 4 for 2-ary, 4-ary, 16-ary respectively
CODE_ARY = 2 ** CODE_PARM
CODE_LEN = 2048
Parity_files_folder = f"Parity_files_{CODE_LEN}"
def main():
    if Parity_files_folder == "Parity_files_8320":
    # code parameter file name
        if CODE_PARM == 1:
            parity_filename = f"{Parity_files_folder}/8320_1280_2ary.dat"
            mapping_filename = "Mapping_files/SignalSet_BPSK-1.txt"
        elif CODE_PARM == 2:
            parity_filename = f"{Parity_files_folder}/4160_640_4ary.dat"
            mapping_filename = "Mapping_files/SignalSet_BPSK-2.txt"
        elif CODE_PARM == 4:
            parity_filename = f"{Parity_files_folder}/2080_320_16ary.dat"
            mapping_filename = "Mapping_files/SignalSet_BPSK-4.txt"
        elif CODE_PARM == 5:
            parity_filename = f"{Parity_files_folder}/1664_256_32ary.dat"
            mapping_filename = "Mapping_files/SignalSet_BPSK-5.txt"
        elif CODE_PARM == 8:
            parity_filename = f"{Parity_files_folder}/1040_160_256aryCode12.dat"
            mapping_filename = "Mapping_files/SignalSet_BPSK-8.txt"
        else:
            raise ValueError(f"Unsupported CODE_ARY: {CODE_ARY}")
    else:
        if CODE_PARM == 1:
            parity_filename = f"{Parity_files_folder}/2048_1024_2ary.dat"
            mapping_filename = "Mapping_files/SignalSet_BPSK-1.txt"
        elif CODE_PARM == 2:
            parity_filename = f"{Parity_files_folder}/1024_512_4aryCode15.dat"
            mapping_filename = "Mapping_files/SignalSet_BPSK-2.txt"
        elif CODE_PARM == 4:
            parity_filename = f"{Parity_files_folder}/512_256_16aryCode17.dat"
            mapping_filename = "Mapping_files/SignalSet_BPSK-4.txt"
        elif CODE_PARM == 8:
            parity_filename = f"{Parity_files_folder}/256_128_256aryCode13.dat"
            mapping_filename = "Mapping_files/SignalSet_BPSK-8.txt"
        else:
            raise ValueError(f"Unsupported CODE_ARY: {CODE_ARY}")
    max_iteration = 50
    k_2 = 320  # inner code length of data bits
    # DNA channel parameters
    # PEs = np.linspace(0.36, 0.3, 4)  # different base error rates
    PEs = [0.304]

    sequencingDepths = [10]  # sequencing depth list
    innerRedundancy = 114  # total redundancy for inner code
    repeat_times = 1000  # repeat simulation times for each sequencing depth

    codec = fftqspa.BCJRQSPA(parity_filename, max_iteration, mapping_filename)

    n_info = codec.info_bits_len()
    n_code = codec.code_bits_len()

    if not os.path.exists('results'):
        os.makedirs('results')

    for sequencingDepth in sequencingDepths:
        runs_filename = f'results/codelen_{CODE_LEN}/FFTQSPA_DNA_Channel_{CODE_ARY}-ary_SequencingDepth{sequencingDepth}_InnerRedundancy{innerRedundancy}_runs.csv'
        results_filename = f'results/codelen_{CODE_LEN}/FFTQSPA_DNA_Channel_{CODE_ARY}-ary_SequencingDepth{sequencingDepth}_InnerRedundancy{innerRedundancy}.csv'

        # 如果 runs_filename 不存在，则创建并写入表头
        if not os.path.exists(runs_filename):
            with open(runs_filename, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Run', 'Pe', 'FER'])
        # 如果 results_filename  存在，则追加写入
        for run_idx in range(1, repeat_times + 1):
            print(f"===== SequencingDepth={sequencingDepth}, Run {run_idx}/{repeat_times} =====")
            for Pe in PEs:
                print(f"Simulating for base error rate Pe={Pe}...")
                # 生成320条长度为n_0的随机信息比特
                rng = np.random.default_rng()
                info_bits = rng.integers(0, 2, size=(n_info, k_2), dtype=np.uint8)

                # 分别对每一条信息比特进行编码
                print(f"{CODE_ARY}-ary LDPC encoding...")
                code_bits = np.empty((n_code, k_2), dtype=np.uint8)
                for i in range(k_2):
                    code_bits[:, i] = codec.encoder4bibo(info_bits[:, i])

                # 使用inner code + DNA存储信道复合信道
                # 输出voting scores 为多数投票得分, 也是就在一个簇中1的占比（约等于P(b=1)），loss sequence的得分为0.5。
                voting_scores = DNAChannel(code_bits, Pe, sequencingDepth, innerRedundancy)  # shape: (n_code, k_2)

                # 将 voting score 以 0.5 为阈值转换为2进制，并与内码输入的 比特流比较计算ber
                hard_decision = (voting_scores >= 0.5).astype(np.uint8)
                bit_errors_voting = np.sum(hard_decision != code_bits)
                ber_voting = bit_errors_voting / code_bits.size
                print(f"Voting score hard-decision BER: {ber_voting}")


                # 计算P(b=0)
                rr_bits_prob = 1 - voting_scores  # P(b=0)

                # Break exact 0.5 ties to avoid symmetric wandering in decoding.
                tie_mask = rr_bits_prob == 0.5
                if np.any(tie_mask):
                    jitter = rng.choice(np.array([-1.0, 1.0]), size=np.count_nonzero(tie_mask)) * 1e-6
                    rr_bits_prob[tie_mask] = 0.5 + jitter

                # 为rr_bits_prob加一个微小值，防止出现0或1的概率，导致LLR无穷大
                epsilon = 0.05
                rr_bits_prob = np.clip(rr_bits_prob, epsilon, 1 - epsilon)

                # 分别对每一条码字进行译码
                print(f"{CODE_ARY}-ary LDPC decoding...")
                decoded_bits = np.empty((n_code, k_2), dtype=np.uint8)
                iters = np.empty(k_2, dtype=int)
                for i in range(k_2):
                    decoded_bits[:, i], iters[i] = codec.decode4bibo(rr_bits_prob[:, i])

                # 分别系统化编码：信息位在码字末尾（前面是校验位）
                sys_start = n_code - n_info
                decoded_bits = decoded_bits[sys_start:, :]  # 信息位

                # Debug
                '''
                # 打印出出错的比特在 LDPC 译码时的帧的位置，以及迭代次数。                
                error_positions = []
                error_iters = []
                # 将decoded_bits转置为 (k_2, n_info)，每行对应一帧的解码结果，方便分析LDPC 译码情况
                info_bits_T = info_bits.T  # shape: (k_2, n_info)
                decoded_bits_T = decoded_bits.T  # shape: (k_2, n_info)
                for i in range(k_2):
                    bit_errors = np.sum(decoded_bits_T[i, :] != info_bits_T[i, :])
                    if bit_errors > 0:
                        error_positions.append(i)
                        error_iters.append(iters[i])
                for i, pos in enumerate(error_positions):
                    print(f"Frame error positions (info bit index): {pos}")
                    print(f"Corresponding decoding iterations: {error_iters[i]}")
                    # 保存对应的帧到磁盘，用于调试，原始信息比特、接收的rr_bits_prob、解码后的比特分别保存到一个csv文件中,文件名用run_idx和pos区分.
                    debug_dir = f'results/codelen_{CODE_LEN}/debug_frames'
                    if not os.path.exists(debug_dir):
                        os.makedirs(debug_dir)
                    frame_prefix = f'run{run_idx}_frame{pos}'
                    # 保存原始信息比特
                    info_bits_filename = os.path.join(debug_dir, f'{frame_prefix}_info_bits.csv')
                    np.savetxt(info_bits_filename, info_bits_T[pos, :], fmt='%d', delimiter=',')
                    # 保存接收的rr_bits_prob
                    rr_bits_prob_filename = os.path.join(debug_dir, f'{frame_prefix}_rr_bits_prob.csv')
                    np.savetxt(rr_bits_prob_filename, rr_bits_prob[:, pos], fmt='%.6f', delimiter=',')
                    # 保存解码后的比特
                    decoded_bits_filename = os.path.join(debug_dir, f'{frame_prefix}_decoded_bits.csv')
                    np.savetxt(decoded_bits_filename, decoded_bits_T[pos, :], fmt='%d', delimiter=',')
                    print(f"error bit index in frame {pos}: ", np.where(decoded_bits_T[pos, :] != info_bits_T[pos, :])[0]+1)
                    print(f"Saved debug files for run {run_idx}, frame {pos} to {debug_dir}")
                    print("===============================================")
                '''
                # 计算误比特率(BER)和帧错误率(FER), 以decoded_bits[0, :]为1帧,及info_bits[0, :]为原始信息
                total_bit_errors = 0
                total_frame_errors = 0
                for i in range(n_info):
                    bit_errors = np.sum(decoded_bits[i, :] != info_bits[i, :])
                    total_bit_errors += bit_errors
                    if bit_errors > 0:
                        total_frame_errors += 1
                BER = total_bit_errors / (n_info * k_2)
                FER = total_frame_errors / n_info
                print(f"Run {run_idx}, Pe={Pe}: BER={BER}, FER={FER}")
                with open(runs_filename, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([run_idx, float(Pe), float(FER)])

        fer_accumulator = {}
        with open(runs_filename, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                pe = float(row['Pe'])
                fer = float(row['FER'])
                if pe not in fer_accumulator:
                    fer_accumulator[pe] = []
                fer_accumulator[pe].append(fer)

        with open(results_filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Pe', 'FER'])
            for pe in sorted(fer_accumulator.keys()):
                final_fer = float(np.mean(fer_accumulator[pe]))
                writer.writerow([pe, final_fer])

        print(f"Saved run-wise FER to {runs_filename}")
        print(f"Saved aggregated FER to {results_filename}")
if __name__ == "__main__":
    main()
