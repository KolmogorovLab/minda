# Minda
###### Note: This tool is under active devlopment.

Minda is a tool for evaluating structural variant (SV) callers that
* standardizes VCF records for compatibility with both germline and somatic SV callers,
* benchmarks against a single VCF input file, or
* benchmarks against an ensemble call set created from multiple VCF input files.

## Installation

Clone the repository and install the dependencies via conda:

```
git clone https://github.com:KolmogorovLab/minda
cd minda
conda env create --name minda --file environment.yml
conda activate minda
./minda.py
```

## Quick Usage

Benchmarking several vcfs against a truth set vcf:

```
./minda.py truthset --base truthset.vcf --vcfs caller_1.vcf caller_2.vcf caller_3.vcf --out_dir minda_out
```

Creating an ensemble from several vcfs and benchmarking against ensemble calls:

```
./minda.py truthset --vcfs caller_1.vcf caller_2.vcf caller_3.vcf --out_dir minda_out
```

## Inputs and Parameters

### Required

#### Truthset

```
--out_dir        path to out directory
--base           path of base VCF
--tsv | --vcfs   tsv file path
                    -OR-
                 vcf file path(s)
```
#### Ensemble
```
--out_dir        path to out directory
--tsv | --vcfs   tsv file path
                    -OR-
                 vcf file path(s)
--min_support |  minimumn number of callers required to support an ensemble call
--conditions        -OR-
                 specific conditions to support a call
```

### Optional
```
--bed            path to bed file for filtering records with BedTool intersect
--filter         filter records by FILTER column; default="['PASS']"
--min_size       filter records by SVLEN in INFO column
--tolerance      maximum allowable bp distance between base and caller breakpoint; default=500
--sample_name    name of sample
--vaf            filter out records below a given VAF treshold
--multimatch     allow more than one record from the same caller VCF to match a single truthset/ensemble record
```
##### VCF Input
Minda standardizes input VCFs by decomposing every SV into start and end records. Records are handled in one of two following ways:
<ol>
    <li>For records having a CHROM:POS pattern in the `ALT` field, the `#CHROM` and `POS` fields are considered the start. Minda then searches for the end record matching the `ALT` field among other records. Alternatively, the `MATEID` from the `INFO` field may be used to find the end record. If no end record is found, the details from the `ALT` field are used to create one.  </li>
    <li>All other records Minda considers start records. The corresponding end records use the start `#CHROM` and `POS` is calculated by adding the start `POS` with absolute value of `SVLEN` or is extracted from the `END` integer in the `INFO` field. </li> 
</ol>
Minda has been tested on VCFs produced by

* Severus
* SAVANA
* nanomonsv
* Sniffles2
* cuteSV
* SVIM
* GRIPSS
* manta
* SvABA.
If you encounter issues with these or other VCF files, please [let us know](https://github.com/KolmogorovLab/minda/issues). 

##### TSV Input
The `--tsv` file has one required column and up three columns. The columns should be as follows:
<ol>
    <li>VCF paths (required)</li>
    <li>caller name</li>
    <li>prefix</li>  
</ol>
If a caller name is not provided, the name listed in the source field of the VCF will be used. If more than one VCF with the same caller name is provided, prefixes disambiguate ID and column names in Minda output files. In the case where prefixes are not provided by the user, Minda automatically assigns a letter prefix in ascending alphabetically order (i.e. A, B, C, etc.).

An example of TSV contents:
```
/path/to/severus_ONT.vcf     Severus     ONT
/path/to/severus_PB.vcf      Severus     PB
/path/to/manta.vcf           manta       ILL
```
##### Specific Conditions
The `--conditions` parameter enables specific user-defined conditions to be met for each ensemble call. Input a list in double quotation marks that contains:

<ol>
    <li>a (nested) list of caller names, each name in single quotation marks with prefixes, if necessary</li>
    <li>an operator in single quoation marks</li>
    <li>a number</li>  
</ol>

For example, from the TSV contents above, to require that an ensemble call be one for which both ONT and PB agree, when using `--tsv` input, specify:
```
"[['ONT_Severus', 'PB_Severus'], '>=', 2]"
```
OR when using `--vcfs` or `--tsv` input:
```
"[[caller_names[:2], '>=', 2]"
```

To combine multiple conditions, add `'&'` or `'|'` between each condition.
For example, to require at least one long-read call and one short-read call to agree, specify for `--tsv` input:
```
"[[['ONT_Severus', 'PB_Severus'], '>=', 1], '&', [['ILL_manta'], '==', 1]]"
```
OR for `--vcfs` or `--tsv` input:
```
"[[caller_names[:2], '>=', 1], '&', [caller_names[2:], '==', 1]]"
```
##### VAF Filtering
###### Note: This requires preprocessing of VCF file.
To run Minda with the `--vaf` parameter, ensure the VCF files have a `VAF` value in the INFO field.  

## Output Files
Both `truthset` and `ensemble` output:
* tp.tsv for each caller
* fp.tsv for each caller
* fn.tsv for each caller
* support.tsv - lists which callers called which truthset/ensemble records
* results.txt - for each caller, lists the overall precision, recall, F1 scores, as well as the number of TP, FN, FP calls overall and by SVTYPE and SVLEN
* removed_records.txt - list of caller IDs of records not evaluated after removing singletons and filtering by FILTER, SVLEN, VAF

`ensemble` also outputs:
* ensemble.vcf

License
-------

Severus is distributed under a BSD license. See the [LICENSE](LICENSE) for details.

Citation
-------
Ayse Keskus, Asher Bryant, Tanveer Ahmad, Byunggil Yoo, Sergey Aganezov, Anton Goretsky, Ataberk Donmez, Lisa A. Lansdon, Isabel Rodriguez, Jimin Park, Yuelin Liu, Xiwen Cui, Joshua Gardner, Brandy McNulty, Samuel Sacco, Jyoti Shetty, Yongmei Zhao, Bao Tran, Giuseppe Narzisi, Adrienne Helland, Daniel E. Cook, Andrew Carroll, Pi-Chuan Chang, Alexey Kolesnikov, Erin K. Molloy, Irina Pushel, Erin Guest, Tomi Pastinen, Kishwar Shafin, Karen H. Miga, Salem Malikic, Chi-Ping Day, Nicolas Robine, Cenk Sahinalp, Michael Dean, Midhat S. Farooqi, Benedict Paten, Mikhail Kolmogorov. "Severus: accurate detection and characterization of somatic structural variation in tumor genomes using long reads." medRxiv 2024, https://doi.org/10.1101/2024.03.22.24304756.

Credits
-------

Minda is being developed in the Kolmogorov Lab at the National Cancer Institute.

Key contributors:

* Asher Bryant
* Ayse Keskus
* Mikhail Kolmogorov

---
### Contact
If you experience any problems or would like to make a suggestion, please submit an [issue](https://github.com/KolmogorovLab/minda/issues).
To contact the developer directly, email asher.bryant@nih.gov.

