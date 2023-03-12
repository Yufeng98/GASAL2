import os
import sys
def get_cells(seqlen, reflen, r):

    # enforce reflen > seqlen
    if seqlen > reflen: 
        reflen, seqlen = seqlen, reflen

    w = r*2 + 1

    # whole matrix contained within band
    if r > reflen:
        return seqlen * reflen
    # elif r > seqlen:
    #     cells = seqlen * reflen - (reflen - r)*(reflen - r)//2
    # else:
    #     cells = seqlen * reflen - (reflen - r)*(reflen - r)//2 - (seqlen - r)*(seqlen - r)//2

    # top left band fills matrix left
    if r > seqlen:
        cells = seqlen * r + (seqlen*(seqlen+1))/2
        if reflen < seqlen + r: # lower right triangle cut off
            x = seqlen + r - reflen
            cells -= (x*(x+1))/2
        return cells

    # top left band leaves matrix left
    cells = w * seqlen - (r*(r+1))/2 # top left triangle cut off
    if reflen < seqlen + r: # lower right triangle cut off
        x = seqlen + r - reflen
        cells -= (r*(r+1))/2 - (x*(x+1))/2
    else:
        cells -= (r*(r+1))/2
    return cells



band_radius = 100

# parse input file
# in_fn = "input/bsw_147_1m_8bit_input.txt"
in_fn = sys.argv[1]
in_fh = open(in_fn, 'r')
lines = [line.rstrip() for line in in_fh.readlines()]
inputs = [tuple(lines[i:i+3]) for i in range(0, len(lines), 3)]
scores, refs, seqs = list(map(list, zip(*inputs)))

# create output files
ref_fn = 'input/bsw_refs.fasta'
seq_fn = 'input/bsw_seqs.fasta'
log_fn = 'input/bsw_log.txt'
if os.path.isfile(ref_fn): os.remove(ref_fn)
if os.path.isfile(seq_fn): os.remove(seq_fn)
if os.path.isfile(log_fn): os.remove(log_fn)
ref_fh = open(ref_fn, 'a')
seq_fh = open(seq_fn, 'a')
log_fh = open(log_fn, 'a')

# reformat output, counting total cells
cells = 0
for i, (ref_str, seq_str) in enumerate(zip(refs, seqs)):

    # int to ACGT
    print(f">ref{i}", file=ref_fh)
    print(f">seq{i}", file=seq_fh)
    try:
        ref_bases = ''.join(["ACGTN"[int(i)] for i in ref_str])
        print(ref_bases, file=ref_fh)
        seq_bases = ''.join(["ACGTN"[int(i)] for i in seq_str])
        print(seq_bases, file=seq_fh)
    except IndexError:
        print("ERROR: Unexpected base.")
        print("reference:", ref_str, "\nsequence:", seq_str)

    # count cells
    new_cells = get_cells(len(seq_bases), len(ref_bases), band_radius)
    cells += new_cells
    print(f'lengths: {len(seq_bases)} '
          f'({band_radius}, {len(ref_bases)}) -> {new_cells}', file=log_fh)

# final output stats for GASAL2 and CUPS calculation
print(f'max ref len: {max(map(len, refs))}', file=log_fh)
print(f'max seq len: {max(map(len, seqs))}', file=log_fh)
print(f'total cells: {cells}', file=log_fh)

# close files
in_fh.close()
ref_fh.close()
seq_fh.close()
log_fh.close()
