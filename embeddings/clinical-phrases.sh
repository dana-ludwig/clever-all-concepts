#!/bin/bash
sed -e "s/’/'/g" -e "s/′/'/g" < /share/pi/stamang/data/mimic/testnotes_formatted.txt | tr -c "A-Za-z-+\/'_ \n" " " > /share/pi/stamang/data/mimic/word2vecOutput/notes-norm0
time word2vec/./word2phrase -train /share/pi/stamang/data/mimic/word2vecOutput/notes-norm0 -output /share/pi/stamang/data/mimic/word2vecOutput/notes-norm0-phrase0 -threshold 200 -debug 2
time word2vec/./word2phrase -train /share/pi/stamang/data/mimic/word2vecOutput/notes-norm0-phrase0 -output /share/pi/stamang/data/mimic/word2vecOutput/notes-norm0-phrase1 -threshold 100 -debug 2
tr A-Z a-z < /share/pi/stamang/data/mimic/word2vecOutput/notes-norm0-phrase1 > /share/pi/stamang/data/mimic/word2vecOutput/notes-norm1-phrase1
time word2vec/./word2vec -train /share/pi/stamang/data/mimic/word2vecOutput/notes-norm1-phrase1 -output /share/pi/stamang/data/mimic/word2vecOutput/vectors-phrase.bin -cbow 1 -size 200 -window 10 -negative 25 -hs 0 -sample 1e-5 -threads 10 -binary 1 -iter 15
word2vec/./distance /share/pi/stamang/data/mimic/word2vecOutput/vectors-phrase.bin
