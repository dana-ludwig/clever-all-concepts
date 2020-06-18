import os, sys, glob
from weakLabelRuleFcns import *

print "Weakly labeling instances and sorting output, by time and for each patient..."

#+++++++++++++++++++++++++++++++++++++++++++++++++++++
# you will need to tweak this for your project
proj = "output" 
criteria = "ALCOHOL-ABUSE"
tag = "ALCOHOL-ABUSE"
projpath = "/share/pi/stamang/people/dana/clever-test/"+proj
fins = glob.glob(projpath+"/ptseq/pt*.txt")
opath = projpath+"/label/"+criteria+"/"

#+++++++++++++++++++++++++++++++++++++++++++++++++++++

fout_pos = open(opath+"allPos"+criteria+".txt","w")
fout_neg = open(opath+"allNonPos"+criteria+".txt","w")
for fin in fins:
  #print fin
  tmp = fin.strip(".txt")
  tmp = tmp.split("pt")
  #print tmp
  tmpf = "pt"+tmp[2]
  #print "PT:",tmpf
  fpos = open(opath+tmpf+"-pos.txt","w")
  fneg = open(opath+tmpf+"-neg.txt","w")
  fcron = open(opath+tmpf+"-cron.txt","w")
  #print "\n"
  with open(fin) as f:
    for line in f:
      sem = 1
      tmp = line.strip()
      label = assignWeakLabel(tmp)
      print label, tmp
      tmpe = tmp.split("|")
      cid = tmpe[0]
      tseq = tmpe[1]
      longseq = tmpe[2]
      tterm = tmpe[3]
      pid = tmpe[4]
      nid = tmpe[5]
      ntype = tmpe[6]
      time = tmpe[7]
      year = tmpe[8]
      tclass = tmpe[9]
      tags = tmpe[14]
      snippet = tmpe[len(tmpe)-1]
      #print label, tmp
      sum_out = label[0]+"|"+label[1]+"|"+pid+"|"+year+"|"+cid+"|"+time+"|"+ntype+"|"+nid+"|"+tterm+"|"+snippet
      long_out =  label[0]+"|"+label[1]+"|"+tmp
      if label[0] == "POSITIVE":
        sem = applyBlackList(snippet,tterm)
        if sem == 1:
          print >> fout_pos, long_out
          print >> fpos, sum_out
          print >> fcron, sum_out
        else: print "FILTERED"
      if label[0] != "POSITIVE" or sem == 0:
        if sem == 0:
          long_out = "NEGATIVE|FILTER|"+tmp
          sum_out = "NEGATIVE|FILTER|"+pid+"|"+year+"|"+cid+"|"+time+"|"+ntype+"|"+nid+"|"+tterm+"|"+snippet
          print "NEW NEG:",long_out
        print >> fout_neg, long_out
        print >> fneg, sum_out
        print >> fcron,sum_out
  fpos.close()
  fcron.close()
  fneg.close()

fout_pos.close()
fout_neg.close()

