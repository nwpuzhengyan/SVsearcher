import pysam
import argparse
from pyfaidx import Fasta
import time

parser = argparse.ArgumentParser()
parser.add_argument("hg19_bam", help='bam file')
args = parser.parse_args()
localtime = time.asctime(time.localtime(time.time()))
print (localtime)
merge_distance=1000

samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
w= open('result2.txt', 'w')

SV_LEN=30
#detect long DEL
samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
long_DEL_record=[]
for r in samfile:
    if r.mapq > 20 and r.is_supplementary == False:
        sign = ""
        for tag in r.tags:
            if str(tag[0])=="SA":
                sign="SA"
                supple_align=tag
        # judge the two alignment
        if sign == "SA":
            # check primary alignment
            chr1 = r.reference_name
            cigar1 = r.cigarstring
            read = r.query_sequence
            other_align = supple_align[1].split(";")
            if r.is_reverse == True:
                direct1='-'
            else:
                direct1 = '+'

            supple=None
            MAPQ_value=0
            i=0
            while i < (len(other_align) - 1):
                m = other_align[i].split(",")
                if m[2]!= direct1 and int(m[4])>MAPQ_value:
                    supple = m
                    MAPQ_value=int(m[4])
                i+=1
            if supple==None:
                continue
            chr2 = supple[0]
            # judge the primary alignment and supplementary alignment pos
            if r.reference_start <= int(supple[1]) and chr1 == chr2:
                pos1 = r.reference_end
                pos2 = int(supple[1])
                cigar2 = supple[3]
                i=0
                middle_len=0
                while i < len(cigar2):
                    num = ''
                    while cigar2[i].isdigit():
                        num = num + cigar2[i]
                        i += 1
                    num = int(num)
                    if cigar2[i] == 'M' or cigar2[i] == 'D':
                        middle_len+=num
                    i+=1
                pos2 = pos2+middle_len

                DEL_len=pos2-pos1
                if DEL_len >= SV_LEN:
                    long_DEL_record.append([chr1, pos1, pos2])
            elif chr1 == chr2:
                pos1 = r.reference_start
                pos2 = int(supple[1])
                cigar2 = supple[3]
                DEL_len = pos1-pos2
                if DEL_len >= SV_LEN:
                    long_DEL_record.append([chr1, pos2, pos1])

res=[]
if len(long_DEL_record):
    intervals = list(sorted(long_DEL_record))
    '''for p in intervals:
        point.write(str(p[0])+':'+str(p[1])+'-'+str(p[2])+'\n')'''
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][2] - intervals[0][1]
    read_num = 1
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][2] - intervals[i][1]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (1.2 * newSV_len) and SV_len >= (0.8 * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
            read_num += 1
        else:
            res.append([chr, low, high, read_num])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][2] - intervals[i][1]
            read_num = 1
    res.append([chr, low, high, read_num])
long_DEL=[]
for site in res:
    if site[3]>=3:
        long_DEL.append(site)
samfile.close()
#merge long DEL and short DEL.
DEL=long_DEL

res=[]
if len(DEL):
    intervals = list(sorted(DEL))
    '''for p in intervals:
        point.write(str(p[0])+':'+str(p[1])+'-'+str(p[2])+'\n')'''
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][2] - intervals[0][1]
    read_num = intervals[0][3]
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][2] - intervals[i][1]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (1.2 * newSV_len) and SV_len >= (0.8 * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
            read_num += intervals[i][3]
        else:
            res.append([chr, low, high, read_num])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][2] - intervals[i][1]
            read_num = intervals[i][3]
    res.append([chr, low, high, read_num])

for site in res:
    if site[3]>20 and (site[2]-site[1])<100000:
       w.write(str(site[0]) + ' ' + str(site[1]) + ' ' + str(site[2])+ '\n')
samfile.close()


#I add modified tep to modify the DEL pos. and add soft clip read step.