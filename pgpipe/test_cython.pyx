

def getAlleleCountArray(char * instr, int slen):
    cdef int i = 0;
    cdef int field_count = 0;
    cdef bint active_field = 0;
    cdef int ref_total = 0;
    cdef int alt_total = 0;
    cdef int miss_total = 0;
    cdef char v = 0
    cdef int d[256]
    for i in range(256): d[i] = 0
    for i in range(slen):
        v = instr[i]
        if v == '\t':
            field_count += 1
            if field_count >= 9:
                active_field = 1
        if v == ';':
            active_field = 0
        if v == '|' or v == '/':
            continue
        if active_field:
            d[v] += 1
            #if v == '0':
            #    ref_total += 1
            #elif v == '1':
            #    alt_total += 1
            #elif v == '.':
            #    miss_total += 1
    return d
    #print (d['0'],d['1'],'0')

def getAlleleCountArraySub(char * instr, int slen, list idxlist, int ilen):
    cdef int i = 0;
    cdef int field_count = 0;
    cdef int iidx = 0;
    cdef bint active_field = 0;
    cdef int d[256];
    cdef char v = 0;
    for i in range(256): d[i] = 0
    for i in range(slen):
        v = instr[i]
        if v == '\t':
            field_count += 1
            if field_count >= 9:
                if idxlist[iidx] == field_count-9:
                    active_field = 1;
                    iidx += 1
                else:
                    active_field = 0
        elif v == ';':
            active_field = 0
        elif v == '|' or v == '/':
            continue
        if active_field:
            d[v] += 1
    return d

         

