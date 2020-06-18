#!/usr/bin/env python
'''
@athors: Suzanne Tamang, Tanya Podchiyska
@date: 06/14/2020
@version: 0.9.2
@status: Development
'''
import glob

POSITIVE="POSITIVE"
NEGATIVE="NEGATIVE"
FILTERED="FILTERED"
DOT="B-DOT"
NEGEX="B-NEGEX"
HYP="B-HYP"
FAM="B-FAM"
HX="B-HX"
SCREEN="B-SCREEN"
ALCOHOL_ABUSE="ALCOHOL-ABUSE"


def gettagterm(tag,dictionary):
  tmp = tag.split(":")
  tterm = dictionary[tmp[1]]
  tmpt = tterm[0].strip()
  if tmpt==":": tterm[0] = "colon"
  elif tmpt==".": tterm[0] = "period"
  elif tmpt=="/": tterm[0] = "slash"
  elif tmpt==",": tterm[0] = "comma" # check ! and ?
  elif tmpt==";": tterm[0] = "semicolon"
  elif tmpt=="!": tterm[0] = "exclamationPoint"
  elif tmpt=="?": tterm[0] = "questionMark"
  taginfo = tterm[0]+":"+tag
  #print tterm[0], tag, taginfo
  return taginfo


def getTagseq(taginfo,window,tclass):
  wmax = int(window)
  wmin = -wmax
  fullseq = ""
  truncseq = ""
  tmptags = taginfo.split("|")
  #print "TAGINFO:",taginfo
  sem = 0
  #print "----------------------"
  for tag in tmptags:
    #print tag
    tmp = tag.split(":")
    tmpclass = tmp[1]
    tmppos = int(tmp[4])
    #print "TPOS:",str(tmppos),tmpclass
    if tmppos < 0:
      if fullseq == "": 
        fullseq = tmpclass
      else: 
        fullseq = fullseq+"_"+tmpclass
      if abs(tmppos) <= wmax:
        #print "good",truncseq
        if truncseq == "":
          truncseq = tmpclass
        else: 
          truncseq = truncseq+"_"+tmpclass
    if tmppos > 0 and sem == 1:
      fullseq = fullseq+"_"+tmpclass
      if tmppos <= wmax:
        truncseq = truncseq+"_"+tmpclass
    if tmppos > 0 and sem == 0:
      if fullseq == "":
        fullseq = "#"+tclass+"#"+"_"+tmpclass
      else:
        fullseq = fullseq+"_"+"#"+tclass+"#"+"_"+tmpclass
      if truncseq == "":
        truncseq = "#"+tclass+"#"+"_"+tmpclass
      else:
        truncseq = truncseq+"_"+"#"+tclass+"#"+"_"+tmpclass
      sem = 1
  #print fullseq, truncseq
  #print tag, taginfo                                                                                                                   
  return [fullseq,truncseq]


def getTerminology(dictname):
  termDict = {}
  print dictname
  with open(dictname) as f:
    for line in f:
      tmp = line.split("|")
      tid = tmp[0].strip()
      term = tmp[1].strip()
      tclass = tmp[2].strip()
      termDict[tid]=[term,tclass]
  return termDict
        

#loads all note metadata for patient subset
def loadSelectNoteMetadata(ptList):
# patient_id|note_id|doc_description|age_at_note_DATE_in_days|note_year        
  print "processing notemetadata for ",str(len(ptList)), " patients"
  fname = "/data3/stride6/tp_annotator_notes.txt"
  noteDict = {}
  fout_notemeta = open("/data3/mbc/notemetadata.txt","w")
  with open(fname) as f_in:
    for line in nonblank_lines(f_in):
      tmp = line.split("|")
      pid = tmp[0].strip()
      nid = tmp[1].strip()
      if pid in ptList:
        print >> fout_notemeta, line.strip()
#       with open("/data3/S6/corpus/notes/"+nid) as onconote:
#         fout = open("/data3/oncoshare/oncocorpus/"+nid,"w")
#         print >> fout, onconote
        #fout.close()
        noteDict[tmp[1]]=[tmp[0],tmp[2],tmp[3],tmp[4]]
  print "Total notes: ",len(noteDict)
  return noteDict

def gdbm_shelve(filename, flag="c"):
  return shelve.Shelf(gdbm.open(filename, flag))

def loadMimicNoteMdata(fname, noteDict): 
# PATIENT_ID|NOTE_ID|MIMIC_ID|TIMESTAMP|NOTE_TYPE
  with open(fname) as f_in:
    next(f_in) # skip header
    for line in nonblank_lines(f_in):
      tmp = line.split("|")
      noteDict[tmp[1]]=[tmp[0],tmp[2],tmp[3],tmp[4]]
  print "Total notes: ",len(noteDict)

def loadRISENoteMdata(fname):
# PATIENT_ID|NOTE_ID|MIMIC_ID|TIMESTAMP|NOTE_TYPE
  with open(fname) as f_in:
    next(f_in) # skip header
    for line in nonblank_lines(f_in):
      tmp = line.split("|")
      noteDict[tmp[1]]=[tmp[0],tmp[2],tmp[3],tmp[4]]
  print "Total notes: ",len(noteDict)

def loadSHCNoteMdata(fname):
# PATIENT_ID|NOTE_ID|MIMIC_ID|TIMESTAMP|NOTE_TYPE
  with open(fname) as f_in:
    next(f_in) # skip header
    for line in nonblank_lines(f_in):
      tmp = line.split("|")
      noteDict[tmp[1]]=[tmp[0],tmp[2],tmp[3],tmp[4]]
  print "Total notes: ",len(noteDict)

#loads all note metadata for STRIDE6
def loadNoteMetadata():
# patient_id|note_id|doc_description|age_at_note_DATE_in_days|note_year
  fname = "/data3/stride6/tp_annotator_notes.txt"
  noteDict = {}
  with open(fname) as f_in:
    for line in nonblank_lines(f_in):
      tmp = line.split("|")
      noteDict[tmp[1]]=[tmp[0],tmp[2],tmp[3],tmp[4]]
  print "Total notes: ",len(noteDict)
  return noteDict

def nonblank_lines(f):
  for l in f:
    line = l.rstrip()
    if line:
      yield line

def hasNumbers(inputString):
  return any(char.isdigit() for char in inputString)

# pt list is a selected group of patients
# set ptList =0 to return all patients
def getPids(ptList,ptkey):
  lc = 0
  mimicMap = {}
  with open(ptkey) as f_in:
    for line in nonblank_lines(f_in):
      if lc == 0: 
        lc += 1
        continue
      line = line.strip()
      tmp = line.split("|")
      #print tmp
      if hasNumbers(tmp) == False: continue
      id_pt = tmp[0].strip()
      id_mimic = tmp[1].strip()
      mimicMap[id_pt]=id_mimic
  f_in.close()
  sample = {}
  if ptList == 0:
    print "Total pts:", len(mimicMap)
    return mimicMap
  with open(ptList) as f_in:
    for line in nonblank_lines(f_in):
      tmp = line.strip()
      id_pt = tmp
      x = mimicMap[id_pt]
      sample[tmp]=x
  print "Total Sampled pts:", len(sample)
  return sample





def assignWeakLabel(truncseq, tclass):
  cseq = formatSeq(truncseq,tclass)
  if cseq == [None,None] : 
    return [POSITIVE,tclass]
  sentsem= checkSentence(cseq[0],cseq[1])
  if sentsem == 1: 
    return [POSITIVE,tclass]
  #print "labeling:", cevent
  label = cleverRule(cseq,tclass)
# print "LABEL: ", label
  #print tmp
  return label 

def formatSeq(seq,tclass):
  lseq = None
  rseq = None
  if "_#"+tclass+"#_" in seq: 
    tmp = seq.split("_#"+tclass+"#_")
    lseq = tmp[0].split("_")
    rseq = tmp[1].split("_")
  elif "_#"+tclass+"#" in seq: 
    tmp = seq.split("_#"+tclass+"#")
    lseq = tmp[0].split("_")
  elif "#"+tclass+"#_" in seq: 
    tmp = seq.split("#"+tclass+"#_")
    rseq = tmp[1].split("_")
# print lseq,rseq
  return [lseq,rseq]

def cleverRule(cseq,tclass):
  #print tclass, cseq 
  supress = [NEGEX,HYP,FAM,HX,SCREEN]
  promote = [ALCOHOL_ABUSE]
  # read in class sequence
  if cseq[0] == None:
    llseq = 0
    pre1 = DOT
  else:
    lseq = cseq[0]
    llseq = len(cseq[0])
    #print lseq
  if cseq[1] == None:
    lrseq = 0
    post1 = DOT
  else:
    rseq = cseq[1]
    lrseq = len(rseq)
    #print rseq
  #print llseq,lrseq
  # see if any of the negation or other nonpositive terms are present
  for tag in supress:
    if llseq > 0: 
      pre1 = lseq[llseq-1]
      #print pre1, tag
      if pre1 == tag: 
        return [NEGATIVE,tag] 
    if lrseq > 0:
      post1 = rseq[0]
      if post1 == tag: 
        return [NEGATIVE,tag]
    if llseq > 2:
      pre2 = lseq[llseq-2]
      if pre2 == tag and pre1 != DOT: 
        return [NEGATIVE,tag]
    if llseq > 3 and tag != NEGEX:
      pre3 = lseq[llseq-3] 
      if pre3 == tag and pre1 != DOT:
        return [NEGATIVE,tag] 
  return [POSITIVE,tclass]

def checkSentence(lseq,rseq):
  if lseq == None and rseq == None: 
    return 1
  elif lseq == None:
    if DOT == rseq[0]: 
      return 1
  elif rseq == None:
    if DOT == lseq[len(lseq)-1]: 
      return 1
  elif DOT == lseq[len(lseq)-1] and DOT == rseq[0]:
      return 1
  else: return 0


def applyBlackList(snippet,tterm):
  s = snippet.lower()
  postTterm = ["no evidence of","any evidence for","any evidence of"]
  for modifier in postTterm:
    if modifier+" "+tterm in s:
      print FILTERED," PREFIX:", modifier+" "+tterm,s
      return 0
  blackList = ["lymph nodes positive","status code:", "did not reveal any disease"]
  for phrase in blackList:
    if phrase in s:
      print FILTERED," PHRASE:", phrase,s
      return 0
    else: return 1