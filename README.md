# cosmic_gaps
Annotation of sequencing gaps with number of variants seen in the COSMIC dataset

## To filter the bed file downloaded from cosmic into referral types:

```
python filter_raw_data.py --cosmic_file <cosmic_bed_file> --config <config_file>
```

## To get the number of cosmic variants in the sequencing gaps after running bedtools intersect:

```
python filter_table.py --sampleId <sampleId> --referral <referral> --gaps_path <gaps_path> --bedfile_path <bedfile_path>
```
