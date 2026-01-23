# MSMuTect
Indel, Allele and Mutation caller, specifically designed to call mutations in microsatellite regions using a pair of sequencing files (normal and tumor samples). It can also be adapted for germline MSI analysis
# Installation
There are 3 possible ways to use MSMuTect: docker container, binary executable, and bash script. The docker container and bash script are around 15% faster than the executable. 

## Binary
NOTE: The 4.1 version binary is not available on github yet. It is available on google drive:
https://drive.google.com/file/d/11a7EwXuZp-dKR6O7e09a6V5UB7Lq_g6K/view?usp=drive_link   
When the bleeding branch moves to main, we will upload the 4.1 binary to github

There is a prebuilt x86_64 linux binary available in releases.  
Note: The binary is slightly slower than the 'Local' option.   
Download the binary from the following link:
[fill in link]

## Docker
To build:   
docker build -t msmutect-docker .   
To run:  
docker run --rm -it msmutect-docker [flags]  
It is necessary to mount the directories of the input and output files so the msmutect docker can access them. 
Here is an example:  
docker run -v /home/avraham/locus_file/:/locus_file -v /home/avraham/data/:/bam_files -v .:/output --rm -it msmutect-docker -l /locus_file/GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed -T /bam_files/tumor.bam -N /bam_files/normal.bam -O /output/docker_test -c 64 -m -A -H



## Bash Script
If on a different platform, or to achieve maximum performance (~15% faster than the binary), do the following:  
git clone https://github.com/MaruvkaLab/MSMuTect_4  
cd MSMuTect_4  
pip3 install -r requirements.txt  
bash build_cython.sh  #optional; improves performance. Otherwise, this is the slowest option
When running, use MSMuTect_4/msmutect.sh everywhere the documentation says 'msmutect'

### Windows Users
MSMuTect cannot run on Windows without WSL because one of MSMuTect's 
dependencies (pysam) is linux-only. While there are replacements for pysam, they are not
considered as mature and stable as pysam. Windows users can use WSL to run
MSMuTect, although this will incur a significant performance hit if the files are
under C:/ (/mnt/c) due to WSL's translation layer

# A Note About Versions
The versions are very confused due to a discontinuous history of development. Practically speaking
you should use either version 4.0 (which is fastest) or version (4.1) which is ~2x slower but 
is much more accurate for longer motif repeats. Both have executable binaries available on the github

# Usage
## Locus File
We strongly recommend using one of the precompiled locus files we have available.
For HG38, the following loci file is best: https://drive.google.com/file/d/1IiWVp1uPz00Daax8z4bP-XI2v_DQPGk7/view?usp=sharing
For HG19, the following loci file is best: https://drive.google.com/file/d/1P9W6VBcpc1bfVXXTiyx9EldK2mbdGiv5/view?usp=sharing
For other genome builds, please email us.
They can be subsampled as desired as long as the order is maintained

If want to use your own locus file, there are a couple of steps you must do:
1. First, the loci file must be sorted properly. This should work on all unix systems:   
sort -t $'\t' -k1,1 -k5n,5 -k4nr,4 -V [original loci file] > [new loci file]    
2. If using version>=4.1, you must also correct the loci file so that the motif is written as it
appears on the forward strand. You can do this with ./tests/scratch/post_process_loci_file.py 
(change the paths as needed)

## Running the Software

msmutect [flags]
  
To fully analyze all loci:  
msmutect -T [tumor_bam.bam] -N [normalbam.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -m -A -H  

To find mutations in the most efficient runtime possible:  
msmutect -T [tumor_bam.bam] -N [normalbam.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -m  

To call indels and alleles for an individual file:  
msmutect -S [sequence_file.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -A -H  

To call indels but not alleles for an individual file:  
msmutect -S [sequence_file.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -H  

To call mutations using results from previously called indels (can only use a single core):
msmutect --from_file -l [loci_file.phobos]  -N [normal.hist.tsv] -T [tumor.hist.tsv]  -O [output_prefix] -m -A 

To see all flags, such as changing parameters, outputting vcf files, etc., run 'msmutect --help'

MSMuTect will create temporary files when running, with names like .msmutect_tmp_file_10242_1721809243.1243694_25529.  
It deletes them at the end. If, for some reason, msmutect is interrupted, these files will not be deleted. They can be safely removed 

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
1. Create a python virtual environment with all of the dependencies in requirements.txt and pyinstaller
2. Run create_executable.sh
