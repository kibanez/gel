The structural variation of repeat elements is a relatively unexplored phenomenon within the field of Repeat expansion disorders (REDs). We present our analysis of short tandem repeats (STRs) using our novel software tool, “Repeat Crawler”. This software supplements output from Expansion Hunter (EH) by taking the sequencing output and recording the presence and length of elements that both interrupt and border the disease repeat within the recorded alleles.

Briefly, Expansion hunter outputs a bam file which contains the reads from a locus of interest. The reads are QC'd and annotated according to the sections of a repetitive region that each read covers. The various components of the repetitive region are marked by 0:flanking region, 1:first component and so on... Hence this can be taken advantage when further analysing reads.

Repeat Crawler accepts a list of repeat component structures (specified in a JSON input file) "crawls" along each read in an EH bamlet file and documents the presence/absence and length of each component in the order that they are arranged in the input JSON file. NOTE the phrase "crawls" is used to convey that the program goes from each individual repeat component to the next as it makes it's way along a read. 
Subsequently it tallies the number of found structures across all reads and finds the 2 most common and reports them as the first and second alleles which were then assigned to the EH allele lengths (matched by adding the total length of the repeat components and comparing to the EH repeat sizes).

Exemplary run of repeat crawler script on a list of Expansion hunter BAMlet files:

```
python find_regex_genotypes.py --bam_files list_of_expansion_hunter_bamlet.txt --output THAP11_RC.tsv --span 0-2 --gene THAP11 \
    --count CAG1 CAA1 CAG2 CAA2 CAG3 CAA3 CAG4 CAA4 CAG5 CAA5 CAG6 CAA6 CAG7 --json THAP11_structure.json \
    --must_contain CAG1
```
