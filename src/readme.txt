Operating instructions for CLEVER on MIMIC III notes as case study:

Modify paths in code as needed to accommodate your setup. 

Step 0A. Preprocess notes to generate a note file in structure note_id<tab>note_text -- i.e. the text of each note appears on one line and is preceded by the note identifier and a tab character -- and a pipe-delimited note metadata file, such as testnotesmdata.txt. 
Assumes presence of a file like testnotes.csv, having following structure: SUBJECT_ID,HADM_ID,ICUSTAY_ID,ELEMID,CHARTTIME,REALTIME,CGID,CORRECTION,CUID,CATEGORY,TITLE,TEXT,EXAM_NAME,PATIENT_INFO

processNotes.py is a python script that transforms the MIMIC notes in testnotes.csv to CLEVER format, e.g.:
/share/pi/stamang/data/mimic/testnotes.csv -> /share/pi/stamang/data/mimic/testnotes_formatted.txt 
and generates a note metadata file of PATIENT_ID|NOTE_ID|MIMIC_ID|TIMESTAMP|NOTE_TYPE

Execution:
# update paths in variables cpath, fin, fout_nmeta, fout_notes in preprocessMimic.py to match your setup.
# update parsing of testnotes file, if deviating from above specified structure.

> cd clever-mimic3/src/ 
> python preprocessMimic.py


Step 0B. Terminology construction: build a word embedding model of the corpus that you can query to perfrom corpus-driven term expansion.  
To build a phrase embedding model, compile word2vec and use the clinical_phrases.sh code to normalize the corpus and learn a neural language model of the corpus.

clinical_phrases.sh is a shell script for normalizing the clinical corpus and training word and phrase embeddings, using a cbow model and the word2vec's source code. It will store its output in directory word2vecOutput, so if you do not have that directory yet:

> mkdir /share/pi/stamang/data/mimic/word2vecOutput

Execution:
> cd step0_embeddings/word2vec
> make
# *malloc.c is a non-standard header file and you may need to update the library to stdlib.c or remove it for compilation

> cd ..
# update paths in clinical-phrases.sh to match your setup and ensure it is an executable
> ./clinical-phrases.sh
# To run the model (based on example directory setup)
> word2vec/distance ../../../../data/mimic/word2vecOutput/vectors-phrase.bin


Place your terminology here: /share/pi/stamang/res/mimic/dicts  

Once terminology is built, start at Step 1.
        
Step 1. Run the CLEVER tagger to generate annotations. 

step1_tagger.py is a Python script for step 1.

Execution (assumes you are in src directory):
# if not specifying the target classes as an argument on the command line below, then script annonates all possible target classes in terminology.
# if argument provided for main-targets, then script annotates only the target classes specified on the command line. 

> python step1_tagger.py --lexicon /share/pi/stamang/res/mimic/dicts/termN2C2.txt --section-headers /share/pi/stamang/res/mimic/headers.txt  ALCOHOL-ABUSE --snippet-length 125 --snippets --notes /share/pi/stamang/data/mimic/testnotes_formatted.txt --workers 2 --output /share/pi/stamang/projects/ctCriteria/ants --left-gram-context 3 --right-gram-context 2 


Step 2. Patient-level event labeling, including extracting sequences of patient-level candidate events. Aggregate and sort patient-level candidate events chronologically using mentions from step 1 and note metadata from preprocessing.

step2_weakLabel.py is a python script for step 2 and depends on functions in step2_patientAntsFcns.py.

step2_weakLabelFcns.py contains functions loadMimicNoteMdata(fname) and loadSeqs(seqFiles,noteDict,mbcDict) that are specific to the metadata file built for MIMIC. Included are example functions for other corpora that can be used to tailor functions for parsing note metadata for new corpora.

# update functions loadMimicNoteMdata and loadSeqs in step2_weakLabelFcns.py to suit your metadata
# verify that you have annotations for your target class(es) from step 1, e.g. in directory projects/ctCriteria/ants

> python step2_weakLabel.py -s MIMIC3 --lexicon /home/tanyapod/clever-n2c2-mimic3-local-final/res/mimic/dicts/termSNOMEDCT.txt --note-metadata /home/tanyapod/clever-n2c2-mimic3-local-final/corpus/testnotemdata-test.txt --input /home/tanyapod/clever-n2c2-mimic3-local-final/proj/ctCriteria/ants --ptseq-output /home/tanyapod/clever-n2c2-mimic3-local-final/proj/ctCriteria/ptseq --label-output /home/tanyapod/clever-n2c2-mimic3-local-final/proj/ctCriteria/labels



 





