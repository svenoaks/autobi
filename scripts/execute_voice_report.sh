#!/bin/sh

for mode in 1
do
    for file in ./datasets/bdc-all/*wav
    do
        java -Xmx4096m -cp AuToBI.jar edu.cuny.qc.speech.AuToBI.VoiceReport  $file $mode #> results_$file_$mode.txt
        exit
    done

done