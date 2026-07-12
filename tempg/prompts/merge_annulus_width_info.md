# A script to merge annulus info

Write a script that takes each text file in the directory:

output/annulus/a_<number>.txt

and concatenate it to the file output/summary/<number>_s_.txt. Do it only
if this concatenation was not done already. The concatenation was already done
if the target file contains the string "annulus width". 

The script should be named "merge_annulus_info" and stored in the scripts/ subdirectory.

