exemplary run of repeat crawler script on a list of Expansion hunter BAMlet files:

```
python find_regex_genotypes.py --bam_files list_of_expansion_hunter_bamlets.txt --output THAP11_RC.tsv --span 0-2 --gene THAP11 \
    --count CAG1 CAA1 CAG2 CAA2 CAG3 CAA3 CAG4 CAA4 CAG5 CAA5 CAG6 CAA6 CAG7 --json THAP11_structure.json \
    --must_contain CAG1
```
