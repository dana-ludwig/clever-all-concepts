import os, sys

def assignWeakLabel(cevent):
  sem = 1
  tmp = cevent.split("|")
  sinfo = tmp[0].split("-")
  cid = sinfo[1]
  #print tmp
  #print sinfo
  tinfo0 = ""
  tinfo1 = ""
  tinfo2 = ""
  #print "cid:", cid
  if len(sinfo)>3:
    tinfo0=sinfo[0]
    tinfo1=cid
    tinfo2=sinfo[2]+"-"+sinfo[3]
    sinfo = [tinfo0,tinfo1,tinfo2]
  #print "sinfo:", sinfo
  tclass = sinfo[2].strip()
  cseq = formatSeq(tmp[1],tclass)
  if cseq == [None,None] : 
    return ["POSITIVE",tclass]
  sentsem= checkSentence(cseq[0],cseq[1])
  if sentsem == 1: 
    return ["POSITIVE",tclass]
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
  pos = "POSITIVE"
  neg = "NEGATIVE"
  supress = ["NEGEX","HYP","FAM","HX","SCREEN"]
  promote = ["ALCOHOL-ABUSE"]
  # read in class sequence
  if cseq[0] == None:
    llseq = 0
    pre1 = "DOT"
  else:
    lseq = cseq[0]
    llseq = len(cseq[0])
    #print lseq
  if cseq[1] == None:
    lrseq = 0
    post1 = "DOT"
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
        return [neg,tag] 
    if lrseq > 0:
      post1 = rseq[0]
      if post1 == tag: 
        return [neg,tag]
    if llseq > 2:
      pre2 = lseq[llseq-2]
      if pre2 == tag and pre1 != "DOT": 
        return [neg,tag]
    if llseq > 3 and tag != "NEGEX":
      pre3 = lseq[llseq-3] 
      if pre3 == tag and pre1 != "DOT":
        return [neg,tag] 
  return [pos,tclass]

def checkSentence(lseq,rseq):
  if lseq == None and rseq == None: 
    return 1
  elif lseq == None:
    if "DOT" == rseq[0]: 
      return 1
  elif rseq == None:
    if "DOT" == lseq[len(lseq)-1]: 
      return 1
  elif "DOT" == lseq[len(lseq)-1] and "DOT" == rseq[0]:
      return 1
  else: return 0


def applyBlackList(snippet,tterm):
  s = snippet.lower()
  postTterm = ["no evidence of","any evidence for","any evidence of"]
  for modifier in postTterm:
    if modifier+" "+tterm in s:
      print "FILTERED PREFIX:", modifier+" "+tterm,s
      return 0
  blackList = ["lymph nodes positive","status code:", "did not reveal any disease"]
  for phrase in blackList:
    if phrase in s:
      print "FILTERED PHRASE:", phrase,s
      return 0
    else: return 1
  
