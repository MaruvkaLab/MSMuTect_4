# MSMuTect
Indel, Allele and Mutation caller, specifically designed to call mutations in microsatellite regions using a pair of sequencing files (normal and tumor samples). It can also be adapted for germline MSI analysis  
*Very important note for existing users*: beginning from version 4.1, your loci file must fit very specific parameters, so many loci files that were usable for previous versions of MSMuTect will not work.  
See the [Locus File](#locus-file) section for more details. 
# Installation
There are 3 possible ways to use MSMuTect: docker container, binary executable, and bash script. 

## Binary
There is a prebuilt x86_64 linux binary available in releases.  
Please use the latest version for the best results.

## Docker
To build:
```bash
docker build -t msmutect-docker .
```
To run:
```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v /path/to/inputs:/inputs:ro -v /path/to/output:/output \
  msmutect-docker [flags]
```

It is necessary to mount the directories of the input and output files so the msmutect docker can access them. 
Here is an example:
```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v /home/avraham/locus_file/:/locus_file:ro \
  -v /home/avraham/data/:/bam_files:ro \
  -v "$PWD":/output \
  msmutect-docker \
  -l /locus_file/GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed \
  -T /bam_files/tumor.bam -N /bam_files/normal.bam \
  -O /output/docker_test -c 64 -m -A -H
```



## Bash Script
The simplest way to start:
```bash
git clone https://github.com/MaruvkaLab/MSMuTect_4
cd MSMuTect_4
pip3 install -r requirements.txt
bash build_cython.sh  #optional; improves performance a few percentage points
```
When running, use MSMuTect_4/msmutect.sh everywhere the documentation says 'msmutect'

### Windows Users
MSMuTect cannot run on Windows without WSL because one of MSMuTect's 
dependencies (pysam) is linux-only. While there are replacements for pysam, they are not
considered as mature and stable as pysam. Windows users can use WSL to run
MSMuTect, although this will incur a significant performance hit if the files are
under C:/ (/mnt/c) due to WSL's translation layer


# Usage
## Locus File
It is almost certain that loci files used with previous versions of MSMuTect will *NOT* work here
We strongly recommend using one of the precompiled locus files we have available. Using a subset of one of the listed loci files will also work seamlessly.
For HG38, the following loci file is best: https://drive.google.com/file/d/1IiWVp1uPz00Daax8z4bP-XI2v_DQPGk7/view?usp=sharing  
For HG19, the following loci file is best: https://drive.google.com/file/d/1P9W6VBcpc1bfVXXTiyx9EldK2mbdGiv5/view?usp=sharing  
The HG38 file has an additional field the HG19 file lacks. This field contains information about the allele distribution
in humans and helps improve MSMuTect's specificity.
For other genome builds, please email us.

We discourage making your own loci file because it will also lack the field the HG19 file lacks. If you are determined to do
so regardless, there are a couple of steps you must do:
1. First, the loci file must be sorted properly. This should work on all unix systems:
   ```bash
   sort -t $'\t' -k1,1 -k5n,5 -k4nr,4 -V [original loci file] > [new loci file]
   ```
2. If using version>=4.1, you must also correct the loci file so that the motif is written as it
appears on the forward strand. You can do this with ./tests/testing_utils/post_process_loci_file.py 
(change the paths as needed)
3. You will still lack the allele distribution field. This will probably hurt the tool's specificity, and MSMuTect will emit a warning about that


## Running the Software

```bash
msmutect [flags]
```
You can also run MSMuTect on CRAM files. To do so, you must also supply the --reference_genome_file flag

To fully analyze all loci:
```bash
msmutect -T [tumor_bam.bam] -N [normalbam.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -m -A -H
```

To find mutations in the most efficient runtime possible:
```bash
msmutect -T [tumor_bam.bam] -N [normalbam.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -m
```

To call indels and alleles for an individual file:
```bash
msmutect -S [sequence_file.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -A -H
```

To call indels but not alleles for an individual file:
```bash
msmutect -S [sequence_file.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -H
```

To call mutations using results from previously called indels (can only use a single core):
```bash
msmutect --from_file -N [normal.hist.tsv] -T [tumor.hist.tsv]  -O [output_prefix] -m -A
```

To see all flags, such as changing parameters, outputting vcf files, etc., run 'msmutect --help'

MSMuTect will create temporary files when running, with names like .msmutect_tmp_file_184260_417968458_1791733248.1329718_1793407325  
It deletes them at the end. If, for some reason, msmutect is interrupted, these files are sometimes not deleted. They can be safely removed 

## Output files
MSMuTect, when run end to end, will generate 3 files
1. A [prefix].full.call_counts.tsv  - this holds a summary of the calls MSMuTect made about the given loci for the input files
2. A [prefix].full.mut.tsv file - this holds the algorithmic calls MSMuTect made about every single locus in the loci file for the input files
3. A [prefix].full.mutated_only.mut.tsv - this holds only the loci that were called mutated in the input file

### Understanding the 'Call' Column
M = Mutation  
NM = Not Mutation  
AN = No Alleles. Either the tumor sample or the normal sample lacks alleles  
RR = Reversion to Reference. The normal sample held an alternative allele, and the tumor had a mutation of the reference allele  
LOH = Loss of Heterozygosity. Tumor has only 1 allele  
FFT = Failed Fisher Test. Passed other tests to be called a mutation, but failed the Fisher's exact test of significance      
INS = Insufficient Support. The normal sample has multiple alleles, but one of them has insufficient support, indicating a noisy locus  
TMA = Too Many Alleles. Normal sample has too many alleles, and hence the locus is too noisy to call    

MSMuTect can generate a vcf file with the results in addition to the regular output (a tsv file).   
The vcf file will include every locus that had an alternate allele in either the tumor or the normal sample.   
However, every locus that is not called as a mutation will be marked as filtered in the filter column of the vcf file, with the resulting call (ex. NM, INS, AN, etc.) as the filter


# Publication and Citation
For orginal paper, see 
YE  Maruvka, Mouw K,  et al, Analysis of somatic microsatellite indels identifies driver events in human tumors
This version is known as version 4.1

# Authors
Avraham Kahan, Dr. Yosef Maruvka, Gaia Frant, and the Maruvka Lab at Technion  
For questions, suggestions, or concerns, open an issue on github or email k.avraham@technion.ac.il or yosi.maruvka@bfe.technion.ac.il

# Changes, Pull Requests, and Compiling the Binary Executable (not relevant for most users)
Users should typically install MSMuTect as described in the "Installation" section. However, if you 
wish to make changes, this is also practicable. The code is straightforward and the entry point is in src/Entry/main.py.
There is a test suite in tests.
When making changes, you don't need to recompile the code every time: simply don't run build_cython.sh and msmutect.sh will use the python files directly.  
If you add any interesting features, we would be grateful if you could share them with us by emailing us or opening a pull request.  
If you did run build_cython.sh and want to revert the .pyx files to python files, run reverse_rename.sh. Please ensure 
you erase the previously generated shared-object libraries, or the python interpreter will use those instead of the .py files
If for some reason you would like to compile a binary executable, this is a little tricky and not recommended.
Specifically, pyinstaller can be picky regarding the exact python and system libraries present on the system

The steps to create the binary executable are as follows:   
1. Create a python virtual environment with all of the dependencies in requirements.txt and pyinstaller. To keep the executable
a manageable size, make the virtual environment as minimal as possible
2. Run build_executable.sh
