#!/usr/bin/env python
'''
step1-tagger.py uses CLEVER's terminology and note header file to generate concept sequences and other annotated textual information for expressing CLEVER rules for automatically labeling events documented in clinical text.
input: file paths to the tagging lexicon with word class mappings, the list of clinical note headers, target classes for event extraction, maximum snippet length, directory path to the clinical corpus, number of workers, output folder and size of n-gram context for n-gram feature generation
output: for target mentions detected using a maximum string length, right truncated partial string matching, CLEVER's output files include right and left n-gram features (context_left.tsv, context_right.tsv), candidate event snippets that can be used for additional processing steps such as SNOMED-CT concept extraction (discover.tsv), and CLEVER's extraction files (extraction.tsv). 
*** it is important to note that only the extraction.tsv file is required to develop a rule based extractor.  Additional textual features are provided in the extraction.tsv file, and other step1-tagger.py output; however, they are inteded to be used in the development of statistical extractors trained on a small portion of development data that is labeled by CLEVER during rule execution
@athors: Suzanne Tamang, Tanya Podchiyska
@date: 06/14/2020
@version: 0.9.2
@status: Development
'''

#from __future__ import division  # needed for NLTK 3.5 use with Python 2
#import nltk, pprint
#from nltk import word_tokenize
#import pdb

import sys
import time
from collections import defaultdict
from multiprocessing import Pool, Queue, Process
import Queue as qmod
from argparse import ArgumentParser
#from os import listdir
#from os.path import isfile, join
import os
import codecs
import re
import traceback
import gc
# from keras.preprocessing.text import text_to_word_sequence # requires Python >= 3.5

#nltk.download()

reload(sys)
sys.setdefaultencoding('utf-8')

# need commas+, otherwise commas in terms lead them to be excluded from annotations! e.g. "influenza,"
END_TOKEN = set([".", " ", "?", "!", "\\", "/", "-", ",", "(", ")", "[", "]", "%", "+", ">", "=", "<", "^", "{", "}", "~"]) # consider adding string.punctuation


def read_headers(f):
  header_list = set()
  try:
    with codecs.open(f, "r", "utf-8") as f:
      for line in f:
        header = line.strip().lower()
        if header:
          header_list.add(header)
  except Exception, e:
    print(e)
    print "Headers file %s failed to load." % f
    print(sys.exc_info())
    print(traceback.print_exc())
  return header_list

def read_dict(f):
  terms = []
  try:
    with codecs.open(f, "r", "utf-8") as f:
      for line in f:
        _id, label, _class = line.strip().split("|")
        # print re.findall("[\w']+", label) # problem w/ term: 2739|++++|C0442725
        word_regex = re.compile("[\w']+", re.IGNORECASE)
        _first_word = word_regex.findall(label)[0] if word_regex.findall(label) else label
        term = Term(_id, label, _class, _first_word)
        terms.append(term)
  except Exception, e:
    print(e)
    print "Term file %s failed to load." % f
    print(sys.exc_info())
    print(traceback.print_exc())
  return terms

class Term:
  def __init__(self, _id, label, _class, _first_word):
    self._id = _id
    self.label = label
    self._class = _class
    self._first_word = _first_word

  def __repr__(self):
    return self.label

class NoteExtraction:
  def __init__(self, note):
    self.note = note
    self.targets = []

  def __repr__(self):
    return "targets %s" % self.targets

  def only_longest_targets(self):
    longest = []
    for target_a in self.targets:
      is_longest = True
      for target_b in self.targets:
        if target_a == target_b:
          continue
        if target_a.is_contained_in(target_b):
          is_longest = False
          break
      if is_longest:
        longest.append(target_a)
    self.targets = sorted(longest, key=lambda x: x.offset)

  def add_target(self, target):
    self.targets.append(target)

  def any_targets(self):
    return len(self.targets) > 0

  def match_headers(self, headers):
    max_target = max(map(lambda x: x.offset, self.targets))
    header_text_space = self.note.text[:max_target]
    header_index = list()
    for header in headers:
      index = [m.start() for m in re.finditer(header+":", header_text_space)]
      for i in index:
        header_index.append((i, header))
    header_index = sorted(header_index, key=lambda x: x[0])
    for target in self.targets:
      for (i, header) in header_index:
        if i > target.offset:
          break
        else:
          target.header = (i, header)

  def dump(self, out_extraction, out_discover, snippets, ngram_contexts): # modify
    sorted_targets = sorted(self.targets, key=lambda x: x.offset)
    lefts, rights = [], []
    for i, target in enumerate(sorted_targets):
      left, right = target.dump(i, out_extraction, out_discover, snippets)
      lefts.append(left)
      rights.append(right)
    return lefts, rights

class MainTargetHit:
  def __init__(self, note, offset, size_context, term, context_terms):
    self.note = note
    self.offset = offset
    self.size_context = size_context
    self.term = term
    self.lsnip, self.rsnip = self.extract_target_snips()
    self.context_terms = context_terms
    self.context_hits=[]
    self.top_offset = self.offset + len(self.term.label)
    self.header = None

  def __repr__(self):
    return "term: %s" % self.term

  def dump(self, i, out_extraction, out_discover, snippets):
    header_str = ""
    if self.header:
      header_str = "%s:%s" % (self.header[1], self.header[0])
    line = [self.note._id, "%s:%s" % (self.term._class, self.term._id),
            str(self.offset), header_str]
    context_sorted = sorted(self.context_hits, key=lambda x: x[2])
    for (offset, term, distance) in context_sorted:
      line.append("%s:%s:%d:%s" % (term._class, term._id, offset,distance))
    if snippets:
      line.append(self.note.text[self.linit:self.rend])
    out_extraction.write("\t".join(line) + "\n")
    left_context = self.note.text[self.linit:self.offset + len(self.term.label)]
    rigth_context = self.note.text[self.offset:self.rend]
    out_discover.write(self.term._class + '\t' + "[[ID=%s:%s:L]]\n"%(self.note._id,str(i)))
    out_discover.write(left_context + "\n")
    out_discover.write(self.term._class + '\t' + "[[ID=%s:%s:R]]\n"%(self.note._id,str(i)))
    out_discover.write(rigth_context+ "\n")
    # return context without the labels
    left = self.note.text[self.linit:self.offset]
    left_w_target = self.term._class + '\t' + left if left else left
    right = self.note.text[self.offset + len(self.term.label):self.rend]
    if not right:
      print(right)
    right_w_target = right + '\t' + self.term._class if right else right
    return left_w_target, right_w_target

  def is_contained_in(self, other):
    return other.offset <= self.offset and other.top_offset >= self.top_offset

  def extract_target_snips(self):
    ltext = len(self.note.text)
    linit = self.offset - self.size_context
    rend = self.offset + len(self.term.label) + self.size_context
    if linit < 0:
        linit = 0
    if rend >= ltext:
       rend = ltext
    self.linit = linit
    self.rend = rend
    l = self.note.text[linit:self.offset]
    r = self.note.text[self.offset + len(self.term.label) : rend]
    return l,r

  def add_context(self, term, hit, side):
    if side == 0: #left
      note_offset = self.offset - (len(self.lsnip) - hit)
    else: #right
      note_offset = self.offset + hit + len(self.term.label)
    distance = note_offset - self.offset
    self.context_hits.append((note_offset, term, distance))

  def only_longest_context(self):
    longest = []
    for context_a in self.context_hits:
      (note_offset, term, distance) = context_a
      top_offset = note_offset + len(term.label)
      is_longest = True
      for context_b in self.context_hits:
        if context_a == context_b:
          continue
        if context_b[0] <= note_offset:
          (note_offset_b, term_b, distance_b) = context_b
          top_offset_b = note_offset_b + len(term_b.label)
          if top_offset_b >= top_offset:
            is_longest = False
            break
      if is_longest:
        longest.append(context_a)
    self.context_hits = sorted(longest, key=lambda x: x[0])

  def extract_context_terms(self):
    for i, snip in enumerate([self.lsnip, self.rsnip]):
      snip_lower = snip.lower()
      for term in self.context_terms:
        offset = 0
        lt = len(term.label)
        ls = len(snip)
        if ls >= lt:
          while True:
            try:
              hit = snip_lower.find(term.label, offset)
            except Exception,e:
              print(e)
              print "Failed finding term %s in snip %s." % (term.label, snip_lower)
              print(sys.exc_info())
              print(traceback.print_exc())
              break
            if hit == -1:
              break
            if hit+lt+1 < ls:
              if snip[hit+lt] not in END_TOKEN:
                break
            self.add_context(term, hit, i)
            # print "added context term %s at hit %d" % (term, hit)
            offset = hit + len(term.label)
    if not args.include_shorter:
      self.only_longest_context()

class Note:
  def __init__(self, _id, text):
    self._id = _id
    self.text = text
    self.text_lower = text.lower()

  def extract(self, main_terms, context_terms, size_context, headers):
    nt = NoteExtraction(self)
    word_regex = re.compile("[\w']+", re.IGNORECASE)
    tokens = set(word_regex.findall(self.text_lower))
    # limit targets and context terms to what is likely to be in note
    main_terms_subset1 = []
    main_terms_subset1 = filter(lambda x: x._first_word in tokens, main_terms)
    tokens = []
    #gc.collect()

    main_terms_subset = []
    for term in main_terms_subset1:
      occurrences = self.text_lower.count(term.label)
      if occurrences > 0:
        main_terms_subset.append(term)
    main_terms_subset1 = []
    
    for term in main_terms_subset:
      offset = 0
      additional_context_terms = []
      try:
        while True:
        #while offset < len(self.text):
          hit = self.text_lower.find(term.label, offset)
          if hit == -1:
            break
          #try:
          #  x=self.text[hit+len(term.label)]
          #except Exception as e:
          #  print(e)
          #  print('note: %d; hit: %d; term: ' + term.label, self._id, hit)
          #  print(sys.exc_info())
          #  print(traceback.print_exc())
          #  break
          if offset == 0 :
            # context is all classes other than target
            additional_context_terms = filter(lambda x: x._class != term._class, main_terms_subset)
            additional_context_terms = additional_context_terms + context_terms
          if hit + len(term.label) < len(self.text):  # check whether term is not at the end of the note
            if self.text[hit + len(term.label)] not in END_TOKEN: # parts of words are not terms
              break
          # tokenization challenge:
          # to avoid tagging terms contained within other words, e.g. ring within during,
          # verify word boundary, TODO: evaluate alternative: split on non-alphanumeric characters
          if hit == 0 or self.text[hit-1] in END_TOKEN:            
            target = MainTargetHit(self, hit, size_context, term, additional_context_terms)
            target.extract_context_terms()
            nt.add_target(target)
          offset = hit + len(term.label)
      except Exception as e:
        print(e)
        print(sys.exc_info())
        print(traceback.print_exc())
    if nt.any_targets():
      nt.match_headers(headers)
      if not args.include_shorter:
        nt.only_longest_targets()
      return nt
    return None

def extract_context(words, start, end):
  words = filter(lambda x: len(x) > 0, words)
  step = 1
  if start > end:
    step = -1
  context = []
  for i in range(start, end, +1):
    if i > -1 and i < len(words):
      context.append(words[i])
    else:
      break
  context = " ".join(context).lower()
  return context

class NGramContext:
  def __init__(self, left_size, right_size):
    self.left_size = left_size
    self.right_size = right_size

  def dump_contexts(self, lefts, rights, fleft, fright):
    #left and right do not include the target label
    if self.left_size:
      for left in lefts:
        if left.find("\t") != -1:
          target_class_l = left.strip().split("\t")[0]
          words = left.strip().split("\t")[1].split(" ")
          words = filter(lambda x: len(x) > 0, words)
          start = len(words)
          lcontext = extract_context(words,
                  start - self.left_size, start)
          if lcontext:
            fleft.write(target_class_l + '\t' + lcontext + "\n")
    if self.right_size:
      for right in rights: 
        if right.find("\t") != -1:
          target_class_r = right.strip().split("\t")[1] if right else ""
          words = right.strip().split("\t")[0].split(" ")
          words = filter(lambda x: len(x) > 0, words)
          rcontext = extract_context(words, 0, self.right_size)
          if rcontext:
            fright.write(target_class_r + '\t' + rcontext + "\n")

  def aggregate(self, output_folder):
    onlyfiles = [join(output_folder, f) for f in listdir(output_folder)\
            if isfile(join(output_folder, f))]
    count_left = defaultdict(lambda : 0)
    count_right = defaultdict(lambda : 0)
    for f in onlyfiles:
      d = None
      if "context-left" in f:
        d = count_left
      elif "context-right" in f:
        d = count_right
      else:
        continue
      with open(f) as fin:
        for line in fin:
          d[line.strip()] += 1
    self.dump_stats(join(output_folder, "context-left-stats.tsv"), count_left)
    self.dump_stats(join(output_folder, "context-right-stats.tsv"), count_right)


  def dump_stats(self, fname, counts):
    total = float(sum(counts.values()))
    accumulative = [(f, c, (float(c)/total) * 100) for (f, c) in counts.items()]
    sorted_acc = sorted(accumulative, key=lambda x: x[1], reverse=True)
    with open(fname, "w") as fout:
      fout.write("\n".join(map(lambda x: "%s\t%d\t%.4f" %
          (x[0], x[1], x[2]), sorted_acc)))


def process_note(line, offset_size, snippets, headers, main_terms, context_terms):
  parts = line.split("\t")
  if len(parts) == 1:
    return
  _id = parts[0]
  text = " ".join(parts[1:])
  note = Note(_id, text)
  note_extraction = None
  try:
    note_extraction = note.extract(main_terms, context_terms, offset_size, headers)
  except Exception,e:
    print(e)
    print "Extraction for note id %s failed." % note._id
    print(sys.exc_info())
    print(traceback.print_exc())      
  finally:
    return note_extraction

class ExitProcess:
  print "pass"
  pass

class Batch:
  def __init__(self, queue, snippet_length, snippets, headers, 
              main_terms, context_terms, output_folder, ngram_contexts):
    self.queue = queue
    self.snippet_length= snippet_length
    self.snippets = snippets
    self.headers = headers
    self.main_terms = main_terms
    self.context_terms = context_terms
    self.output_folder = output_folder
    self.ngram_contexts = ngram_contexts

    if isinstance(self.queue, basestring):
      self.notes_file = open(self.queue, "r")

  def next_batch(self):
    if not isinstance(self.queue, basestring):
      try:
        batch = queue.get(True)
        if isinstance(batch, ExitProcess):
          return None
        return batch
      except qmod.Empty:
        print("qmod empty")
        return []
    else:
      if self.notes_file == None:
        return None
      batch = []
      for line in self.notes_file:
        batch.append(line.strip())
      self.notes_file = None
      return batch

  def process(self):
    if isinstance(self.queue, basestring):
      pid = 0
    else:
      pid = os.getpid()
    output_file = codecs.open(
        os.path.join(self.output_folder, "extraction-%d.tsv"%pid), 
        "w", encoding='utf8')
    discover_file = codecs.open(
        os.path.join(self.output_folder, "discover-%d.tsv"%pid), 
        "w", encoding='utf8')
    fcontext_left = codecs.open(
        os.path.join(self.output_folder, "context-left-%d.tsv"%pid), 
        "w", encoding='utf8')
    fcontext_right = codecs.open(
        os.path.join(self.output_folder, "context-right-%d.tsv"%pid), 
        "w", encoding='utf8')

    while True:
      try:
        batch = self.next_batch()
        if batch is None:
          return
        for line in batch:
          try: 
            ext = process_note(line, self.snippet_length,
                              self.snippets, headers,
                              self.main_terms, self.context_terms)
            if ext:
              lefts, rights = ext.dump(output_file, discover_file, 
                                        self.snippets, self.ngram_contexts)
              if self.ngram_contexts:
                ngram_contexts.dump_contexts(lefts, rights,
                                              fcontext_left, fcontext_right)
          except Exception,e:
            print(e)
            print(sys.exc_info())
            print(traceback.print_exc())
            continue 
      except Exception,e:
        print(e)
        print(sys.exc_info())
        print(traceback.print_exc())
        #traceback.print_stack()
        continue

if __name__ == "__main__":
  parser = ArgumentParser()
  parser.add_argument("-o", "--output", dest="output_folder", default=None,
                      help="output folder", metavar="FILE")
  parser.add_argument("-n", "--notes", dest="notes_file", default=None,
                      help="Notes file", metavar="FILE")
  parser.add_argument("-l", "--lexicon", dest="lexicon", default="mbc-dic.txt",
                      help="read word classes from FILE", metavar="FILE")
  parser.add_argument("-w", "--workers", dest="workers", default=2,
                      metavar="N")
  parser.add_argument("-s", "--section-headers", dest="section_headers",
                      default="headers.txt", help="read headers from FILE",
                      metavar="FILE")
  parser.add_argument("-t", "--main-targets", dest="main_targets",
                      action="append", default=[],
                      help=("the word classes to use as a main target "
                            "(can be used multiple times)"),
                      metavar="TARGET") # optional
  parser.add_argument("-ln", "--snippet-length",
          dest="snippet_length", type=int, default=150)
  parser.add_argument('--snippets', dest='snippets', action='store_true')
  parser.add_argument('--shorter-too', dest='include_shorter', action='store_true',
          default=False)
  parser.add_argument('--no-snippets', dest='snippets',  action='store_false')
  parser.add_argument('--left-gram-context', dest='left_gram',default=3)
  parser.add_argument('--right-gram-context', dest='right_gram',default=2)
  args = parser.parse_args()
  args.workers = int(args.workers)

  if not args.output_folder:
    print("Output folder must be provided with -o/--output")
    sys.exit(-1)
  if not os.path.exists(args.output_folder):
    os.mkdir(args.output_folder)
    print("Output folder '%s' created for saving clean data" % (args.output_folder))
  elif not os.path.isdir(args.output_folder):
    print("Please provide a folder with -o/--output")
    sys.exit(-1)
  elif os.listdir(args.output_folder):
    print(os.listdir(args.output_folder))
    print("Please provide an empty folder for storing annotations of target classes. '%s' already contains data."%(args.output_folder))
    sys.exit(-1)

  terms = read_dict(args.lexicon)
  headers = read_headers(args.section_headers)

  if not args.snippets and (args.right_gram > 0 or args.left_gram > 0):
    sys.stderr.write(("If snippets are disabled context "
                      "ngrams cannot be extracted"))
    sys.exit(-1)
  ngram_contexts = None
  
  if args.left_gram:
    args.left_gram = int(args.left_gram)
  if args.right_gram:
    args.right_gram = int(args.right_gram)

  if args.snippets and (args.right_gram > 0 or args.left_gram > 0):
    ngram_contexts = NGramContext(args.left_gram, args.right_gram)

  # alternatively, you can initialize main_targets_index from a file or on the command line
  main_targets_index = set([x._class for x in terms if not(x._class.startswith('B-'))]) 
  #print(" ".join(main_targets_index))
  if len(args.main_targets) > 0:
    main_targets_index = set(map(lambda x: x.strip(),
                                  args.main_targets[0].split(",")))

  main_terms = filter(lambda x: x._class in main_targets_index,terms)
  context_terms = filter(lambda x: x._class not in main_targets_index,terms)
  terms = [] # no longer needed

  if len(main_terms) == 0:
    sys.stderr.write("Main targets not found - exiting")
    sys.exit(-1)
  
  if args.workers > 0:
    queue = Queue(args.workers)
    batch = Batch(queue, args.snippet_length, args.snippets,
                headers, main_terms, context_terms, args.output_folder,
                ngram_contexts)
    pool = Pool(args.workers, batch.process)

    batch = []
    try:
      with open(args.notes_file,"r") as file_notes:  # verify whether rb is not needed
        for line in file_notes:
          try:
            if len(batch) == 100: #5000
              queue.put(batch, True, None)
              batch = []
            batch.append(line.strip())
          except Exception, e:
            print(e)
            print "Batch at length %d." % len(batch)
            print "Failed to load into batch line %s." % line              
            print(sys.exc_info())
            print(traceback.print_exc())
            continue
        if batch:
          queue.put(batch, True, None)

    except Exception,e:
      print(e)
      print "opening notes file %s." % args.notes_file
      print(sys.exc_info())
      print(traceback.print_exc())
      
    time.sleep(5)
    for x in range(args.workers):
      queue.put(ExitProcess())
    queue.close()
  else:
    batch = Batch(args.notes_file, args.snippet_length, args.snippets,
                headers, main_terms, context_terms, args.output_folder,
                ngram_contexts)
    batch.process()
      
  #if ngram_contexts:
  #  ngram_contexts.aggregate(args.output_folder)
