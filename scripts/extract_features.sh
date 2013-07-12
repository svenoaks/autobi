#!/bin/sh


java -cp AuToBI.jar edu.cuny.qc.speech.AuToBI.ExtractAcousticsSS datasets/dirndl/dlf-nachrichten-200703260000.wav 50 300 5 69.41 70.66 75.47 76.03 > yojojoi.txt
./extract_acoustics-ss.pl /Users/vitisoto/Documents/summer/AndrewRosenberg-AuToBI-67ba0b5/datasets/dirndl/dlf-nachrichten-200703260000.wav 50 300 5 69.41 70.66 75.47 76.03
