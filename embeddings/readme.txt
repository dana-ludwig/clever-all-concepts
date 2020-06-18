Instructions for Step 1: Terminology Construction

1) process all notes so that each note appears on one line and is preceded by a tab character and the note identifier.

processNotes.py will transform the MIMIC notes to CLEVER format:
/share/pi/stamang/data/mimic/testnotes.txt -> /share/pi/stamang/data/mimic/textnotes_formatted.txt

> cd clever-n2c2-mimiciii/src/ 
> python processMimic.py

2) build a word embedding model of the corpus that you can query to perfrom corpus-driven term expansion.  This is descriped as Step 1 in the CELVER framework, Terminology Construction.

> cd clever/src/embeddings/word2vec/
> make

*malloc.c is a non-standard header file and you may need to update the library to stdlib.c or remove it for compilation

To build a phrase ebedding model, compile word2vec and use the clinical_phrases.sh code to normalize the corpus and learn a neural language model of the corpus.

clinical_phrases.sh will store its output in directory word2vecOutput:

> mkdir /share/pi/stamang/data/mimic/word2vecOutput

To run the model,

> cd clever/src/embeddings/word2vec
> ./distance /share/pi/stamang/data/mimic/word2vecOutput/vector-phrase.bin

3) Place your terminology here: /share/pi/stamang/res/mimic/dicts

