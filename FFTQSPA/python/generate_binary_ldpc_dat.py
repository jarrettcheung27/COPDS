import argparse
import os


def parse_degree3_file(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    if len(lines) < 3:
        raise ValueError("Input file is too short.")

    # Header
    # Line 0: RowNumber****ColNum******Rank
    # Line 1: row col rank
    header_parts = lines[1].split()
    if len(header_parts) < 3:
        raise ValueError("Header line does not contain row/col/rank.")
    row_num = int(header_parts[0])
    col_num = int(header_parts[1])
    rank = int(header_parts[2])

    # Find the start of row definitions
    try:
        start_idx = lines.index("1's_location_per_row") + 1
    except ValueError as exc:
        raise ValueError("Missing '1's_location_per_row' marker.") from exc

    rows = []
    for i in range(start_idx, len(lines)):
        parts = lines[i].split()
        if len(parts) < 2:
            continue
        row_id = int(parts[0])
        weight = int(parts[1])
        cols = list(map(int, parts[2:2 + weight]))
        if len(cols) != weight:
            raise ValueError(f"Row {row_id} has inconsistent weight data.")
        rows.append((row_id, cols))

    if len(rows) != row_num:
        raise ValueError(f"Expected {row_num} rows but found {len(rows)} rows.")

    # Build matrix
    matrix = [[0] * col_num for _ in range(row_num)]
    for row_id, cols in rows:
        for c in cols:
            matrix[row_id][c] = 1

    return row_num, col_num, rank, matrix


def sys_matrix_gf2(matrix, row_num, col_num):
    # Mimic Qary_Gauss_E::SysMatrix behavior for GF(2)
    tempP = list(range(col_num))
    enc_matrix = [row[:] for row in matrix]
    tempH = [row[:] for row in matrix]

    m_codechk = 0
    for i in range(row_num):
        flag = False
        ii = jj = None
        for j in range(i, col_num):
            for r in range(i, row_num):
                if enc_matrix[r][j] != 0:
                    flag = True
                    ii = r
                    jj = j
                    break
            if flag:
                break

        if not flag:
            break

        m_codechk += 1

        # swap rows i and ii
        if ii != i:
            enc_matrix[i], enc_matrix[ii] = enc_matrix[ii], enc_matrix[i]

        # swap columns i and jj
        if jj != i:
            tempP[i], tempP[jj] = tempP[jj], tempP[i]
            for r in range(row_num):
                enc_matrix[r][i], enc_matrix[r][jj] = enc_matrix[r][jj], enc_matrix[r][i]

        # elimination (GF(2))
        for r in range(row_num):
            if r != i and enc_matrix[r][i] != 0:
                # row_r = row_r + row_i (mod 2)
                row_r = enc_matrix[r]
                row_i = enc_matrix[i]
                for c in range(col_num):
                    row_r[c] ^= row_i[c]

        # normalize row i (already 1)

    # Reorder original matrix columns according to tempP
    permuted_matrix = [[0] * col_num for _ in range(row_num)]
    for r in range(row_num):
        for c in range(col_num):
            permuted_matrix[r][c] = tempH[r][tempP[c]]

    return m_codechk, permuted_matrix, enc_matrix


def matrix_to_sparse_rows(matrix):
    rows = []
    for r, row in enumerate(matrix):
        cols = [c for c, v in enumerate(row) if v != 0]
        rows.append((r, cols))
    return rows


def write_dat(path, row_num, col_num, rank, q_ary, parity_rows, enc_rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write("row_number\n")
        f.write(f"{row_num}\n")
        f.write("col_number\n")
        f.write(f"{col_num}\n")
        f.write("rank\n")
        f.write(f"{rank}\n")
        f.write("ary\n")
        f.write(f"{q_ary}\n")
        f.write("*****1's_location/per_row****\n")
        for row_id, cols in parity_rows:
            f.write(f"{row_id}   {len(cols)}   ")
            for c in cols:
                f.write(f"{c}   1     ")
            f.write("\n")
        f.write("*****1's_location/per_row****\n")
        for row_id, cols in enc_rows:
            f.write(f"{row_id}   {len(cols)}   ")
            for c in cols:
                f.write(f"{c}   1     ")
            f.write("\n")


def write_signalset(path):
    with open(path, "w", encoding="utf-8") as f:
        f.write("***MappingChart***\n")
        f.write("number_of_signal\n")
        f.write("2\n")
        f.write("length_of_signal\n")
        f.write("1\n")
        f.write("signal_set\n")
        f.write("-1\n")
        f.write("1\n")


def main():
    parser = argparse.ArgumentParser(description="Generate binary LDPC .dat and SignalSet_BPSK.txt")
    parser.add_argument("--input", required=True, help="Path to 8320_1280_degree3.txt")
    parser.add_argument("--dat-output", required=True, help="Output .dat file path")
    parser.add_argument("--signal-output", required=True, help="Output SignalSet_BPSK.txt path")
    args = parser.parse_args()

    row_num, col_num, rank, matrix = parse_degree3_file(args.input)
    m_codechk, permuted_matrix, enc_matrix = sys_matrix_gf2(matrix, row_num, col_num)

    # Use computed rank if different
    rank_out = m_codechk if m_codechk > 0 else rank

    parity_rows = matrix_to_sparse_rows(permuted_matrix)
    enc_rows = matrix_to_sparse_rows(enc_matrix[:rank_out])

    write_dat(args.dat_output, row_num, col_num, rank_out, 2, parity_rows, enc_rows)
    write_signalset(args.signal_output)


if __name__ == "__main__":
    main()
