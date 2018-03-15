import sys



sample_count = 20
missing_data_count = [0 for i in range(sample_count+1)]

indel_count = 0

missing_data_per_indiv = [0 for i in range(sample_count)]

prev_miss_pos = 0
prev_full_pos = 0
prev_pos = 0
prev_chrom = ''

in_section = False

for line in sys.stdin:
    if line[0] == '#':
        continue
    la = line.strip().split()
    missing_for_site = False
    for i in range(sample_count):
        b_i = 9+i
        geno = la[b_i]
        if '.' in geno:
            missing_for_site = True
            break
    cur_pos = int(la[1])
    if prev_pos == 0 or prev_chrom != la[0]:
        prev_miss_pos = cur_pos
        prev_full_pos = cur_pos
        prev_chrom = la[0]
        if not missing_for_site:
            in_section = True
        prev_pos = cur_pos
        continue
    if missing_for_site:
        if in_section:
            in_section = False
            diff = prev_pos - prev_full_pos
            if diff > 1000:
                print (la[0],prev_miss_pos,prev_full_pos,prev_pos,cur_pos,diff)
    else:
        if not in_section:
            prev_miss_pos = prev_pos
            prev_full_pos = cur_pos
            in_section = True
    prev_pos = cur_pos
