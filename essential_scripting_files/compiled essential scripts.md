## Codes for extraction information from fasta/bam file
- **To get the read's mate unmapped from bam file** 
~~~
samtools view <file>.bam | cut –f7 | grep –c ‘*’
