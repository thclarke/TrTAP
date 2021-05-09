#!/bin/env python3
import sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[1], "rb")
diff = int(sys.argv[2])
header = samfile.header

list2 = {}
lens = {}
list_gd = {}
tot_1 = {}
tot_2 = {}
tot_3 = {}
minV = {}
tot = []
for i in samfile.fetch(until_eof=True):
    if i.is_paired & i.is_read1:
      if i.query_name in list2:
          if i.has_tag('AS'):
             if minV[i.query_name]:
              if i.get_tag('AS') == minV[i.query_name]:
                list2[i.query_name].append(i.reference_name)
              else:
                if i.get_tag('AS') > minV[i.query_name]:
                  list2[i.query_name] = [i.reference_name]
                  minV[i.query_name] = i.get_tag('AS')
          else:
             list2[i.query_name] = [i.reference_name]
             minV[i.query_name] = i.get_tag('AS')
      else:
         if i.has_tag('AS'):
           minV[i.query_name] = i.get_tag('AS')
           list2[i.query_name] = [i.reference_name]
      if len(i.tags) > 2:
         if i.tags[2][1] < diff:
            if i.query_name in list_gd:
                list_gd[i.query_name].append(i.reference_name)
            else:
                list_gd[i.query_name] = [i.reference_name] 
print("Seq ID\tNorm Count\tUnambig Count\tGood Count\tNorm FPKM\tUnambig FPKM\tGood FPKM")
tot.append(0.0)
tot.append(0.0)
tot.append(0.0)
tot.append(0.0)
for i in list2.keys():
    q = list2[i]
    #print(i)
    #print(q[0])
    if len(q) == 1:
        tot[1] = tot[1] + 1
        tot[2] = tot[2] + 1
        if q[0] in tot_1:
            tot_1[q[0]] = tot_1[q[0]] + 1
        else:
            tot_1[q[0]] = 1
        if q[0] in tot_2:
            tot_2[q[0]] = tot_2[q[0]] +1
        else:
            tot_2[q[0]] = 1
        if i in list_gd:
            tot[3] = tot[3] + 1
            if q[0] in tot_3:
                tot_3[q[0]] = tot_3[q[0]] +1
            else:
                tot_3[q[0]] =1
    else:
        for j in q:
            tot[1] = tot[1] + 1.0/len(q)
            if j in tot_1:
                tot_1[j]= tot_1[j] +1.0/len(q)
            else:
                tot_1[j] = 1.0/len(q)
tot[1] = tot[1] / 1000000000
tot[2] = tot[2] / 1000000000
tot[3] = tot[3] / 1000000000

for i in samfile.references:
    q = 0.0;
    if i in tot_1:
         q = tot_1[i]
    #print(i)
    n1 = 0
    len1 = float(samfile.get_reference_length(i) * 1.0) / 1000.0
    #print(str(len1) + "\t" + str(samfile.lengths[i]))
    if i in tot_3:
        n1 = tot_3[i]
    n3_fpkm = n1 / (tot[3] * len1)
    n1_fpkm = q / (tot[1] * len1)
    if i in tot_2:
        if (q > -1):
            print(str(i)+ "\t"+ str(q)+ "\t"+ str(tot_2[i])+ "\t"+ str(n1) + "\t" + str(n1_fpkm) + "\t" + str((tot_2[i]/(tot[2] * len1))) + "\t" + str(n3_fpkm))
        else:
            print("No_Hit\t"+ str(q)+ "\t"+ str(tot_2[i])+ "\t"+ str(n1) + "\t0\t0\t0")
    else:
        if (q > -1):
            print(str(i)+ "\t"+ str(q)+ "\t0\t"+ str(n1) + "\t" + str(n1_fpkm) + "\t0\t" + str(n3_fpkm))
        else:
            print("No_Hit"+ "\t"+ str(q)+ "\t0\t"+ str(n1) +"\t0\t0\t0")
