# MEI-VCF-to-FASTA
This tool converts the INFO column from a MELT vcf into fasta format. The script also requires either the single **transposon fasta file** used in MELT analysis with option **-f** or the designation of one of the hg19 transposon families supplied by MELTv2.14 (**-m Alu**).

# Usage<br/>
```
usage: Convert_MELT_vcf_to_fasta.py [-h] -i INPUT VCF [-o OUTPUT FASTA]
                                    [-m MELT HG19 TRANSPOSON]
                                    [-f TRANSPOSON FASTA]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT FASTA, --out OUTPUT FASTA
                        output FASTA file. Default will have the same prefix as vcf
  -m MELT HG19 TRANSPOSON
                        MELT HG19 human transposon files (Alu, LINE1, or SVA)
  -f TRANSPOSON FASTA, --fasta TRANSPOSON FASTA
                        transposon fasta used for MELT analysis

required named arguments:
  -i INPUT VCF, --in INPUT VCF
                        name/path to VCF

```

Please email me if you have any questions: jfeusier@genetics.utah.edu
