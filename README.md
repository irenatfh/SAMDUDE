# SAMDUDE
SAMDUDE is a genomic sequence denoiser that operates on aligned SAM files. It both denoises individual bases in reads as well as updates the corresponding quality scores.

Irena Fischer-Hwang ihwang@stanford.edu

## Installation

SAMDUDE requires no additional software beyond Python 3+ and NumPy.

## Usage

The input SAM files must be sorted. Denoising affects only the read and quality score strings; all headers, identifiers and other fields (flags, CIGAR strings, etc.) will be kept intact.

To denoise a sam file, use the following command:

```
python3 run_denoiser.py noisy_file.sam denoised_file.sam
```
