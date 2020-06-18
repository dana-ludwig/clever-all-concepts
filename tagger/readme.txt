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

In the projects directory, create output directories for each of the target classes.  For example, the pqrs anntotation output directories are located here: /share/pi/stamang/projects/ctCriteria

The example below targets the "ALCOHOL-ABUSE" class and writes annotations to the /share/pi/stamang/projects/ctCriteria/ALCOHOL-ABUSE output directory for further processing.

Execution (assumes you are in src directory):
# if you do not already have a directory ants/ for target class annotations, then create it:
> mkdir /share/pi/stamang/projects/ctCriteria/ants

# if you have already run the tagger for this target class, rename or remove the annotations directory with the target class name, so that the tagger can run again, i.e. for this example rename or remove ALCOHOL-ABUSE directory in /share/pi/stamang/projects/ctCriteria/ants.

# if not specifying the target classes as an argument on the command line below, then set them in variable main_targets_index

++++++++++++++++++++++++++++++++++
python2 tagger.py --lexicon /share/pi/stamang/cie/clever-demo/res/termN2C2.txt --section-headers /share/pi/stamang/cie/clever-demo/res/headers.txt --main-targets ALCOHOL-ABUSE --snippet-length 125 --snippets --notes /share/pi/stamang/cie/clever-demo/notes/testnotes_formatted.txt --workers 4 --output /share/pi/stamang/cie/clever-demo/etoh/ants/ALCOHOL-ABUSE --left-gram-context 3 --right-gram-context 2
++++++++++++++++++++++++++++++++++

Step 2. Extract sequences of patient-level candidate events. Aggregate and sort patient-level candidate events for target classes chronologically using mentions from step 1 and note metadata from preprocessing.

step2_patientAnts.py is a python script for step 2 and depends on functions in step2_patientAntsFcns.py.

step2_patientAntsFcns.py contains functions loadMimicNoteMdata(fname) and loadSeqs(seqFiles,noteDict,mbcDict) that are specific to the metadata file built for MIMIC. Included are example functions for other corpora that can be used to tailor functions for parsing note metadata for new corpora.

Execution:
# if you don't already have a ptseq directory, create it:
> mkdir /share/pi/stamang/projects/ctCriteria/ptseq

# update paths for variables dictfile, project, ppath, noteMdata and set target_class(es) in step2_patientAnts.py to match your setup. The target classes are the same ones used in Step 1. ptKey can be used to limit patients to a preidentified set. Without it, the script will run over any patient with annotations for the specified target classes.

# update functions loadMimicNoteMdata and loadSeqs in step2_patientAntsFcns.py to suit your metadata
# verify that you have annotations for your target class(es) from step 1 in directory projects/ctCriteria/ants

> nohup python step2_patientAnts.py &


Step 3. Patient-level event labeling.

Execution:
# if you don't already have a labels directory, create it:
> mkdir /share/pi/stamang/projects/ctCriteria/labels

# update variables proj, criteria, tag, projpath, fins, opath, fout_pos, fout_neg in script step3_weakLabelALCOHOL-ABUSE.py or customize step3_weakLabelTARGET-CLASS.py. If you modified code to accommodate your metadata structure, edit these scripts accordingly.
# these scripts have a dependency on step3_weakLabelRuleFcns.py 

# if you don't already have a directory for storing the labels of each target class, create it:
> mkdir /share/pi/stamang/projects/ctCriteria/labels/ALCOHOL-ABUSE

> nohup python step3_weakLabelALCOHOL-ABUSE.py &


 





