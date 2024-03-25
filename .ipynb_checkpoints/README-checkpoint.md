# Minda
Minda is tool for benchmarking structural variant (SV) callers that:
* requires no prepocessing of VCF files
* standardizes VCF records for compatibility with both germline and somatic SV callers
* benchmarks against a single VCF input file or ensemble call set (created from multiple VCF input files)

## Installation

Clone the repository and install the dependencies via conda:

```
git clone https://github.com:KolmogorovLab/minda
cd Severus
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
##### TSV Input
The `--tsv` file has one required column and up three columns. The columns should be as follows:
<ol>
    <li>VCF paths (required)</li>
    <li>caller name</li>
    <li>prefix</li>  
</ol>
If a caller name is not provided, the name listed in the source field of the VCF will be used. If more than one VCF with the same caller name is provided, prefixes disambiguate ID and column names in Minda output files. In the case where prefixes are not provided by the user, Minda automatically assigns a letter prefix in ascending alphabetically order (i.e. A, B, C, etc).

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

Using the the TSV contents above, for example, to require that an ensemble call be one for which both ONT and PB agree, specify for `--tsv` input:
```
"[['ONT_Severus', 'PB_Severus'], '>=', 2]"
```
OR for `--vcfs` or `--tsv` input:
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

## Output Files
Both `truthset` and `ensemble` output:
* TP TSV for each caller
* FP TSV for each caller
* FN TSV for each caller
* Support TSV - lists which callers called which truthset/ensemble records
* Results - for each caller, lists the overall precision, recall, F1 scores, as well as the number of TP, FN, FP calls overall and by SVTYPE and SVLEN
* Removed records - list of caller IDs of records not evaluated after removing singletons (breakends not having a mate) and filtering by FILTER, SVLEN, VAF

`ensemble` also outputs:
* Ensemble VCF

License
-------

Severus is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

Minda is being developed in the Kolmogorov Lab at the National Cancer Institute.

Key contributors:

* Asher Bryant
* Ayse Keskus
* Mikhail Kolmogorov

---
### Contact
For suggestions, reporting bugs or help, please submit an [issue](https://github.com/KolmogorovLab/minda/issues).
To contact the developer, email asher.bryant@nih.gov.

