#!/usr/bin/env python

'''
@athors: Suzanne Tamang, Tanya Podchiyska
@date: 06/14/2020
@version: 0.9.2
@status: Development
'''

import os, sys, operator
from argparse import ArgumentParser
from operator import itemgetter
from step2_weakLabelRuleFcns import *

MIMIC = 'MIMIC3'
RISE = 'RISE'
SHC = 'SHC'

POSITIVE="POSITIVE"
NEGATIVE="NEGATIVE"
FILTERED="FILTERED"

if __name__ == "__main__":
  parser = ArgumentParser() 
  parser.add_argument("-lo", "--label-output", dest="label_output_folder", default=None,
                        help="weak label directory", metavar="DIR") 
  parser.add_argument("-i", "--input", dest="ants_folder", default=None,
                        help="Extraction directory", metavar="DIR")
  parser.add_argument("-m", "--note-metadata", dest="noteMdata", default="testnotemdata.txt",
                        help="read note metadata from FILE", metavar="FILE")    
  parser.add_argument("-l", "--lexicon", dest="lexicon", default="termSNOMEDCT.txt",
                        help="read word classes from FILE", metavar="FILE")
  parser.add_argument("-t", "--target-class", dest="target_class",
                        action="append", default=[],
                        help=("the word classes to use as a main target "
                              "(can be used multiple times)"),
                        metavar="TARGET")
  parser.add_argument("-s", "--data-source", dest="data_source", default="MIMIC3",
                        help=("data source used, such as MIMIC3, RISE, or SHC"),
                        metavar="SOURCE")
  # limiting patients to a preidentified set of say patient_id|mimic_id :
  parser.add_argument("-pt", "--ptKey", dest="ptKey", default="mimicKey.txt",
                        help="read preidentified patient list from FILE", metavar="FILE") 

  args = parser.parse_args()
  # args.workers = int(args.workers) # reminder to add workers and multiprocessing

  if not args.ants_folder:
    print("Extraction folder must be provided with -i/--input")
    sys.exit(-1)
  elif not os.listdir(args.ants_folder):
    print("Extraction folder '%s' is empty. Please provide a folder that contains extractions."%(args.ants_folder))
    sys.exit(-1)

  if not args.lexicon:
    print("Terminology must be provided with -l/--lexicon")
    sys.exit(-1)
    
  if not args.noteMdata:
    print("Note metadata folder must be provided with -m/--note-metadata")
    sys.exit(-1)
    
  if not args.data_source:
    print("Data source must be provided with -s/--data-source and a corresponding parsing function must have been built.")
    sys.exit(-1)
    
  if not args.label_output_folder:
    print("Weak label output folder must be provided with -lo/--label-output")
    sys.exit(-1)
  if not os.path.exists(args.label_output_folder):
    os.mkdir(args.label_output_folder)
    print("Weak label output folder '%s' created for saving clean data"%(args.label_output_folder))
  elif not os.path.isdir(args.label_output_folder):
    print("Please provide a folder with -lo/--label-output")
    sys.exit(-1)
  elif os.listdir(args.label_output_folder):
    print("Please provide an empty folder for storing weak labels.'%s' already contains data."%(args.label_output_folder))
    sys.exit(-1)

  #m3pts = getPids(0,ptKey)                    #pt ids for subset
  termDict = getTerminology(args.lexicon)          #load terminology

  noteDict = {}
  if args.data_source == MIMIC:
    loadMimicNoteMdata(args.noteMdata, noteDict)    #load note metada w/format specific to corpus                         
  #s6:patient_id|note_id|doc_description|age_at_note_DATE_in_days|note_year
  #M3: PATIENT_ID|NOTE_ID|MIMIC_ID|TIMESTAMP|NOTE_TYPE
  elif args.data_source == RISE:
    loadRISENoteMdata(args.noteMdata, noteDict)
  elif args.data_source == SHC:
    loadSHCNoteMdata(args.noteMdata, noteDict)

  # alternatively, you can initialize main_targets_index from a file or on the command line
  target_class = set(map(itemgetter(1), filter(lambda x: not(x[1].startswith('B-')), termDict.values())))
  #print target_class
  if len(args.target_class) > 0:
    target_class = set(map(lambda x: x.strip(), args.target_class[0].split(",")))

  if len(target_class) == 0:
    sys.stderr.write("Target classes not found - exiting")
    sys.exit(-1) 

  antsFiles = args.ants_folder + "/extraction*.tsv"
  
  for name in glob.glob(antsFiles):
    print "Loading: ",name
    fin = open(name,"r") # utf-8?
    for line in fin:
      line = line.strip()
      tmp = line.split("\t")
      cols = len(tmp)
      nid = tmp[0].strip()
      if nid not in noteDict.keys(): 
        continue
      mdata = noteDict[nid]
      tmpt = tmp[1].split(":")
      tclass = tmpt[0]
      ttid = tmpt[1]
      tpos = tmp[2]
      
      tmpterm = termDict[ttid]
      tterm = tmpterm[0]

      tags = []
      tmp_tags = []
      tag_offests = []
      sem = 0
     
      pid = mdata[1]    # patient ID in the corpus, i.e. MIMIC id
      sid = mdata[0]    # assigned pid
      doc_desc = mdata[3]
      tage = mdata[2]   # timestamp in MIMIC; age in other cases, etc.
      age = str(tage)
      nyear = "XXXX"            
      snippet = tmp[cols-1].strip()
      #print "snippet of note_id", nid,  "for tclass", tclass, ": ", snippet
      
      if cols < 4: 
        continue 
      if sem != 1:
        # section header processing
        tmpTestHead = tmp[3].split(":")   
        if tmpTestHead[0].islower() and tmpTestHead != "": 
          htmp = tmp[3].split(":")
          head = htmp[0]
          hpos = htmp[1]
          if cols>5: tmptags = tmp[4:cols-1] 
          else: tmptags = []
        else:
          head = "UK"
          hpos = "NULL"
          tmptags = tmp[4:cols-1]
      tmp_key = "S"+"-"+pid+"-"+tclass
      if tmptags == []:
        tmpstr="NONE"
        tagseqs = ["NONE","NONE"]
      else:
        tagterm = gettagterm(tmptags[0],termDict)
        tmpstr = tagterm
        if len(tmptags) > 1:
          for i in range(1,len(tmptags)):
            tmp_item = tmptags[i]
            tagterm = gettagterm(tmp_item,termDict)
            tmpstr = tmpstr+"|"+tagterm
        tagseqs = getTagseq(tmpstr,"125",tclass)
      truncseq = tagseqs[1]
      fullseq = tagseqs[0]
      sinfo = "|".join([tmp_key,truncseq,fullseq,tterm,sid,nid,doc_desc,age,nyear,tclass,ttid,tpos,head,hpos,tmpstr,"SNIPPET: "+snippet])
      #print tagseqs
      #print sinfo
      #print tmptags
      #print tmp
      
      label = assignWeakLabel(truncseq,tclass) # e.g. label = ["POSITIVE",tclass]
      long_out =  "".join(["[",label[0],",",label[1],"]|",sinfo])
      if label[0] == POSITIVE:
        sem = applyBlackList(snippet,tterm)
        if sem != 1:
          long_out="FILTERED " + long_out
      if label[0] != POSITIVE or sem == 0:
        if sem == 0:
          long_out = "NEGATIVE FILTER " + long_out

      if long_out.startswith("[POSITIVE"):
        #fout_pos = open(args.label_output_folder+"/"+tclass+"-pos.txt","a+")
        fout_pos = open(args.label_output_folder+"/allPos.txt","a+")
        print >> fout_pos, long_out     
        fout_pos.close()
      else: 
        #fout_neg = open(args.label_output_folder+"/"+tclass+"-nonPos.txt","a+")
        fout_neg = open(args.label_output_folder+"/allNonPos.txt","a+")
        print >> fout_neg, long_out
        fout_neg.close()
  print "fini!"
