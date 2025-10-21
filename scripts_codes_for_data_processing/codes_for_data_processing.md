Codes for routine tasks related to sequncing data
- *To extract genes information for certain regions of genome based on "start" and "end" column/subsetting genes list from a file for certain genomic region.**
~~~
## Below codes are for running as a script
#!/bin/bash
var=$@
echo $var
var2=$(($var+5000))
var1=$(($var-5000))
echovar "$var1;$var2";
##### based on two conditions in two columns
cat /home/farid/Vflay_genes_list | awk -F "\t" -v start="$var1" -v end="$var2" '{ if(($3>=start) && ($4<=end)) { print } }' 
##### three conditions in two columns
awk -F " " -v start="$var1" -v end="$var2" '{ if(($1=="SOVchr2") && (($3>=start && $3<=end) && ($4>=end))) { print } }' 
condition c
awk -F " " -v start="$var1" -v end="$var2" '{ if(($1=="SOVchr2") && (($3<=start) && ($4>=start && $4<=end))) { print } }' 
~~~
- *To get the read's mate unmapped from bam file** 
~~~
samtools view <file>.bam | cut –f7 | grep –c ‘*’
~~~
