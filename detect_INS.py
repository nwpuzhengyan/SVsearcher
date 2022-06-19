import pysam
import argparse
from pyfaidx import Fasta
import time
import copy

parser = argparse.ArgumentParser()
parser.add_argument("hg19_bam", help='bam file')
args = parser.parse_args()
localtime = time.asctime(time.localtime(time.time()))
print (localtime)
merge_distance=1000
low_rate=0.8
high_rate=1.2
mini_SV_LEN=40

samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
w= open('result2.txt', 'w')
#point= open('site2.txt', 'w')
record=[]
#detect the candidate DEL
for r in samfile.fetch():
  if r.mapq>20 and r.is_supplementary==False:
      chr = r.reference_name
      start_pos=r.reference_start
      cigar=r.cigarstring
      read_str=r.query_sequence
      i=0
      pos=start_pos
      pos1=0
      pos2=0
      read_i = 0
      while i<len(cigar):
          num=''
          while cigar[i].isdigit():
              num=num+cigar[i]
              i+=1
          num=int(num)
          if cigar[i]=='M':
              pos=pos+num
              read_i = read_i + num
          elif cigar[i]=='I':
              pos1 = pos
              pos2 = pos + num
              if num >= 30:
                  INS_seq = read_str[read_i:read_i + num]
                  record.append([chr, pos1, pos2, INS_seq])
              read_i = read_i + num
          elif cigar[i]=='D':
              pos=pos+num
          elif cigar[i]=='S':
              read_i = read_i+num
          i+=1
samfile.close()
#merge the candidate DEL
localtime = time.asctime( time.localtime(time.time()) )
print (localtime)
res=[]

if len(record):
    intervals = list(sorted(record))
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = len(intervals[0][3])
    SV_sequence = intervals[0][3]
    read_num = 1
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][2] - intervals[i][1]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr:
            low = min(intervals[i][1],low)
            high = max(intervals[i][2],high)
            SV_sequence = intervals[i][3]
            read_num += 1
        else:
            res.append([chr, low, high, read_num, SV_sequence])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = len(intervals[i][3])
            SV_sequence = intervals[i][3]
            read_num = 1
    res.append([chr, low, high, read_num, SV_sequence])
localtime = time.asctime( time.localtime(time.time()))
print (localtime)
newres=[]
for site in res:
    if int(site[3])>=3 and len(site[4])>=30:
        newres.append([site[0],site[1],site[2],site[4]])
res_region=newres
#modified the DEL pos(For same read, there are two alignments.)

newDEL=[]
samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
for region in res_region:
    chr = region[0]
    sv_pos1 = region[1]
    sv_pos2 = region[2]
    print(chr, sv_pos1)
    read_del=[]
    zero_DEL_num = 0
    coverage=0
    for r in samfile.fetch(chr, sv_pos1, sv_pos2+1):
        read_tag=r.tags
        align_num=1
        for tag in r.tags:
            if str(tag[0])=="SA":
                sign="SA"
                supple_align=tag
                other_align = supple_align[1].split(";")
                other_align = other_align[:-1]
                align_num=align_num+len(other_align)

        if r.mapq > 20 and r.is_supplementary == False and \
                r.reference_start<(sv_pos1 - 1000) and r.reference_end>(sv_pos2 + 1000):
            coverage+=1
            range1=sv_pos1 - 1000
            range2 = sv_pos2 + 1000
            DEL_num=0
            DEL_len=0
            DEL_len_10bp = 0
            DEL_pos=[]
            start_pos = r.reference_start
            cigar = r.cigarstring
            read_str = r.query_sequence
            read_name=r.query_name
            i = 0
            pos = start_pos
            pos1 = 0
            pos2 = 0
            read_i = 0
            while i < len(cigar):
                num = ''
                while cigar[i].isdigit():
                    num = num + cigar[i]
                    i += 1
                num = int(num)
                if cigar[i] == 'M':
                    pos = pos + num
                    read_i = read_i + num
                elif cigar[i] == 'I':
                    pos1 = pos
                    pos2 = pos+num
                    INS_seq = read_str[read_i:read_i + num]
                    if pos1 > range1 and pos1 < range2:
                        if num >= 40:
                            DEL_pos.append([chr, pos1, pos1,num])
                            DEL_num += 1
                            DEL_len += num
                    read_i = read_i + num
                elif cigar[i] == 'D':
                    pos = pos + num
                elif cigar[i] == 'S':
                    read_i = read_i + num
                i += 1
            if DEL_num!=0:
                read_del.append([DEL_num, DEL_len, DEL_pos])
            else:
                zero_DEL_num+=1
            if coverage>150:
                break

    if len(read_del)==0 or float(len(read_del))/(float(len(read_del))+float(zero_DEL_num))<0.2:
        continue
    read_del = sorted(read_del, key=lambda x: x[1])
    noise_len = len(read_del) / 20
    if noise_len != 0:
        read_del = read_del[noise_len:-noise_len]

    sv_len_num = []
    sv_len = read_del[0][1]
    read = read_del[0]
    num = 1
    i = 1
    while i < (len(read_del)):
        if read_del[i][1] == sv_len:
            if read_del[i][0] < read[0]:
                read = read_del[i]
            num += 1
        else:
            sv_len_num.append([sv_len, num, read])
            sv_len = read_del[i][1]
            read = read_del[i]
            num = 1
        i += 1
    sv_len_num.append([sv_len, num, read])

    cengci_cluster = []
    if len(sv_len_num) <= 2:
        cengci_cluster.append(sv_len_num)
    while len(sv_len_num) > 2:
        distance = []
        distance_rate = []
        i = 1
        while i < len(sv_len_num):
            node_dis = sv_len_num[i][0] - sv_len_num[i - 1][0]
            node_dis_rate = float(sv_len_num[i][0] - sv_len_num[i - 1][0]) / float(sv_len_num[i - 1][0])
            distance.append(node_dis)
            distance_rate.append(node_dis_rate)
            i += 1

        i = 1
        min_dis = distance[0]
        min_dis_rate = distance_rate[0]
        sign = 0
        while i < len(distance):
            if distance[i] < min_dis:
                min_dis = distance[i]
                min_dis_rate = distance_rate[i]
                sign = i
            i += 1

        if min_dis_rate<=0.2 or min_dis<=200:
            node_support_num = sv_len_num[sign][1] + sv_len_num[sign + 1][1]
            node_len = (float(sv_len_num[sign][0]) * sv_len_num[sign][1] + float(sv_len_num[sign + 1][0]) *
                        sv_len_num[sign + 1][1])
            node_len = node_len / node_support_num
            if sv_len_num[sign][2][0] < sv_len_num[sign + 1][2][0]:
                node_read = sv_len_num[sign][2]
            elif sv_len_num[sign][2][0] > sv_len_num[sign + 1][2][0]:
                node_read = sv_len_num[sign + 1][2]
            else:
                if abs(sv_len_num[sign][2][1] - node_len) < abs(sv_len_num[sign + 1][2][1] - node_len):
                    node_read = sv_len_num[sign][2]
                else:
                    node_read = sv_len_num[sign + 1][2]
            del sv_len_num[sign:sign + 2]
            sv_len_num.insert(sign, [node_len, node_support_num, node_read])
            sv_len_num_2 = copy.deepcopy(sv_len_num)
            cengci_cluster.append(sv_len_num_2)
        else:
            break

    if len(sv_len_num) == 1:
        chr = sv_len_num[0][2][2][0][0]
        del_pos1 = sv_len_num[0][2][2][0][1]
        del_pos2 = sv_len_num[0][2][2][-1][2]
        del_len = sv_len_num[0][0]
        newDEL.append([chr, del_pos1, del_pos2, sv_len_num[0][1], del_len])
    elif len(sv_len_num) == 2:
        cengci_cluster.pop()

        sv_len_num = sorted(sv_len_num, key=lambda x: x[1])
        sv_len_num.reverse()
        sv_len_num = sv_len_num[0:2]

        total_num = sv_len_num[0][1] + sv_len_num[1][1]
        support_read1=sv_len_num[0][1]
        support_read2=sv_len_num[1][1]
        rate1 = float(sv_len_num[0][1]) / float(total_num)
        rate2 = float(sv_len_num[1][1]) / float(total_num)
        del1=[sv_len_num[0][2][2][0][0],sv_len_num[0][2][2][0][1],sv_len_num[0][2][2][-1][2],sv_len_num[0][1],sv_len_num[0][0]]
        del2 = [sv_len_num[1][2][2][0][0], sv_len_num[1][2][2][0][1], sv_len_num[1][2][2][-1][2],sv_len_num[1][1],sv_len_num[1][0]]
        average_len = sv_len_num[0][1] * sv_len_num[0][0] + sv_len_num[1][1] * sv_len_num[1][0]
        average_len = float(average_len) / float(total_num)
        len_dif = abs(sv_len_num[1][0] - sv_len_num[0][0])
        while rate1<0.2 or rate2<0.2:
            if len(cengci_cluster)==0:
                rate1 = 0.5
                rate2 = 0.5
                len_dif = 0
                break
            cluster=sorted(cengci_cluster[-1], key=lambda x: x[1])
            total_num=cluster[-1][1]+cluster[-2][1]
            support_read1 = cluster[-1][1]
            support_read2 = cluster[-2][1]
            rate1=float(cluster[-1][1])/float(total_num)
            rate2 = float(cluster[-2][1]) / float(total_num)
            del1 = [cluster[-1][2][2][0][0], cluster[-1][2][2][0][1], cluster[-1][2][2][-1][2],cluster[-1][1],
                    cluster[-1][0]]
            del2 = [cluster[-2][2][2][0][0], cluster[-2][2][2][0][1], cluster[-2][2][2][-1][2], cluster[-2][1],
                    cluster[-2][0]]
            average_len = cluster[-1][1] * cluster[-1][0] + cluster[-1][1] * cluster[-1][0]
            average_len = float(average_len) / float(total_num)
            len_dif = abs(cluster[-1][0] - cluster[-2][0])
            cengci_cluster.pop()
            if cengci_cluster==0:
                break

        len_dif_value = 20 + float(0.005 * average_len)
        if len_dif > len_dif_value:
            if rate1 >= 0.2 and  support_read1 >= 3:
                newDEL.append(del1)
            if rate2 >= 0.2 and  support_read2 >= 3:
                newDEL.append(del2)
        else:
            if (rate1+rate2)>=0.2 and (support_read1+support_read2)>=3:
                newDEL.append(del1)

samfile.close()
print (localtime)
short_DEL=newDEL

mini_len=300
#detect long DEL
samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
long_DEL_record=[]
for r in samfile.fetch():
    if r.mapq > 20 and r.is_supplementary == False:
        read_name=r.query_name
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
            supple = other_align[0].split(",")
            i = 1
            while i < (len(other_align) - 1):
                m = other_align[i].split(",")
                if int(m[4]) > int(supple[4]):
                    supple = m
                i += 1
            chr2 = supple[0]
            # judge the primary alignment and supplementary alignment pos
            if r.reference_start <= int(supple[1]) and chr1 == chr2:
                pos1 = r.reference_end
                pos2 = int(supple[1])
                direct2 = supple[2]
                cigar2 = supple[3]
                s2 = cigar2.split("S")[0]
                if 'M' in s2:
                    s2 = 0
                else:
                    s2 = int(s2)
                if r.is_reverse == True:
                    direct1 = '-'
                else:
                    direct1 = '+'
                read_pos1 = r.query_alignment_end
                read_pos2 = s2
                INS_len = read_pos2 - read_pos1 - (pos2 - pos1)
                INS_seq = read[read_pos1:read_pos1+INS_len]
                if INS_len >= mini_len and direct2 == direct1:
                    long_DEL_record.append([chr1,min(pos1,pos2),max(pos1,pos2),INS_len])
            elif chr1 == chr2:
                pos1 = r.reference_start
                pos2 = int(supple[1])
                direct2 = supple[2]
                cigar2 = supple[3]
                r2_len = 0
                i = 0
                while i < len(cigar2):
                    num = ''
                    while cigar2[i].isdigit():
                        num = num + cigar2[i]
                        i += 1
                    num = int(num)
                    if cigar2[i] == 'M' or cigar2[i] == 'D':
                        r2_len = r2_len + num
                    i += 1
                s2 = 0
                if cigar2[-1] == 'S' or cigar2[-1] == 'H':
                    s2 = num
                pos2 = pos2 + r2_len
                if r.is_reverse == True:
                    direct1 = '-'
                else:
                    direct1 = '+'
                read_pos2 = len(read) - s2
                read_pos1 = r.query_alignment_start
                INS_len = read_pos1 - read_pos2-(pos1-pos2)
                INS_seq = read[read_pos2:read_pos2+INS_len]
                if INS_len >= mini_len and direct2 == direct1:
                    long_DEL_record.append([chr1,min(pos1,pos2),max(pos1,pos2),INS_len])

res=[]

if len(long_DEL_record):
    intervals = list(sorted(long_DEL_record))
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][3]
    #SV_sequence = intervals[0][4]
    read_num = 1
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][3]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (high_rate * newSV_len) and SV_len >= (low_rate * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][3]
            #SV_sequence = intervals[i][4]
            read_num += 1
        else:
            res.append([chr, low, high, read_num, SV_len])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][3]
            #SV_sequence = intervals[i][4]
            read_num = 1
    res.append([chr, low, high, read_num, SV_len])
long_DEL=[]
for site in res:
    if site[3]>=3:
        long_DEL.append(site)
samfile.close()
#merge long DEL and short DEL.
DEL=short_DEL+long_DEL

res=[]
if len(DEL):
    intervals = list(sorted(DEL))
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][4]
    #SV_sequence = intervals[0][5]
    read_num = 1
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][4]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (high_rate * newSV_len) and SV_len >= (low_rate * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][4]
            #SV_sequence = intervals[i][5]
            read_num += 1
        else:
            res.append([chr, low, high, SV_len])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][4]
            #SV_sequence = intervals[i][5]
            read_num = 1
    res.append([chr, low, high, SV_len])
samfile.close()

candidate_res=[]
for site in res:
    if site[3]>=40:
        candidate_res.append([str(site[0]),int(site[1]),int(site[2]),int(site[3]), 0])

total_DEL=candidate_res

if len(total_DEL)!=0:
    w.write(str(total_DEL[0][0]) + ' ' + str(total_DEL[0][1]) + ' ' + str(total_DEL[0][2]) + ' ' + \
                 str(total_DEL[0][3]) + ' ' + str(total_DEL[0][4]) + '\n')
i=1
while i<len(total_DEL):
    if total_DEL[i-1][0]==total_DEL[i][0] and (int(total_DEL[i][1])-int(total_DEL[i-1][2]))<1000 and\
            int(total_DEL[i-1][3])<10000 and int(total_DEL[i][3])<10000:
        dif_rate = float(abs(total_DEL[i - 1][3] - total_DEL[i][3])) / float(min(total_DEL[i - 1][3], total_DEL[i][3]))
        if dif_rate < 0.2 and total_DEL[i - 1][3] > 200 and total_DEL[i][3] > 200:
            total_DEL[i - 1][4] = 2
            total_DEL[i][4] = 2
        else:
            total_DEL[i - 1][4] = 1
            total_DEL[i][4] = 1
            w.write(str(total_DEL[i][0]) + ' ' + str(total_DEL[i][1]) + ' ' + str(total_DEL[i][2]) + ' ' + \
                     str(total_DEL[i][3]) + ' ' + str(total_DEL[i][4]) + '\n')
    else:
        w.write(str(total_DEL[i][0]) + ' ' + str(total_DEL[i][1]) + ' ' + str(total_DEL[i][2]) + ' ' + \
                 str(total_DEL[i][3]) + ' ' + str(total_DEL[i][4]) + '\n')
    i+=1



#I add modified tep to modify the DEL pos. and add soft clip read step.