# Chimpanzee

1.0.0

A pipeline for te-gene chimera pair detection from Nanopore Direct cDNA Sequencing or Direct RNA Sequencing data.

## Usage

Clone this repository into your desired directory on HARDAC:

```bash
git clone https://github.com/subunit/chimpanzee.git
```

Edit the parameters in the `chimpanzee.sh` script:

```bash
samplename=         # sample name
experimenttype=     # DCS or DRS
species=            # fly or mouse (currently not used in pipeline)

inputdir=           # path to directory with input fastq files
referencegenome=    # path to reference genome
referencete=        # path to te genome
geneanno=           # path to annotation file in bed12 format
refseqid=           # path to refseq-symbol conversion file

minimap2=           # path to minimap2
runpychopper=       # Y or N (N recommended)

netid=              # Duke NetID
minimumhit=         # Minimum hits of te-gene chimera to count as significant
```

Make scripts executable:

```bash
chmod u+x chimpanzee.sh
chmod u+x chimpanzee_align.sh
```

Run pipeline by submitting `chimpanzee.sh` as a wrapper script:

```bash
sbatch --mem=64G --cpus-per-task=4 --ntasks=1 --wrap="./chimpanzee.sh"
```
