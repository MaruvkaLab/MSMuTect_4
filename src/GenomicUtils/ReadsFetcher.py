# cython: language_level=3
import pysam
from typing import List, Optional
from pysam import AlignmentFile, AlignedSegment
from pysam.libcalignmentfile import IteratorRowRegion

from src.GenomicUtils.AlignmentFlags import FLAG_OPTIONS


def open_alignment_file(path: str, reference: Optional[str] = None) -> AlignmentFile:
    """
    Open a BAM or CRAM file with pysam.

    CRAM is reference-based, so the reference genome is required to decode its
    records and is passed through to pysam. BAM is self-contained, so the
    reference is ignored for .bam files.
    """
    if path.endswith(".cram"):
        return pysam.AlignmentFile(path, reference_filename=reference)
    elif path.endswith(".bam"):
        # BAM is self-contained; the reference is ignored.
        return pysam.AlignmentFile(path)
    else: # should never happen
        raise RuntimeError(f"Unrecognized file type: {path}. This error should not happen. Please report this error to the developers.")



class ReadsFetcher:
    """
    ReadsFetcher fetches reads from the BAM file
    ReadsFetcher stores what it last returned to a query, and checks these results for candidates for the current query
    It can miss reads if the loci file is improperly sorted. These parts of the program are fairly tightly coupled, unfortunately
    """
    def __init__(self, alignment_file: str, start_chromosome: str, reference_genome_file: Optional[str] = None):
        self.BAM_handle = open_alignment_file(alignment_file, reference=reference_genome_file)
        self.chromosome = start_chromosome
        self.chromosome_prefix = self.get_prefix()
        self.last_extracted_reads = []
        self.reads_iterator: IteratorRowRegion = self.BAM_handle.fetch(f"{self.chromosome_prefix}{self.chromosome}", start=10_000, multiple_iterators=False)
        self.last_unmapped_read: AlignedSegment = self.get_next_mapped_read()  #  will hold last non mapped read to previous query so will not be lost in the iterator

    def get_prefix(self) -> str:
        """ returns the prefix of contigs in the BAM (ex. Chr, nothing [''], or chr)"""
        prefixes = ['', 'chr', 'Chr']
        for prefix in prefixes:
            try:
                _ = self.BAM_handle.fetch(f"{prefix}{self.chromosome}", start=10_000, multiple_iterators=False)
                return prefix
            except ValueError:  # different prefix
                continue
        print("UNKNOWN CHROMOSOME PREFIX FOUND: PLEASE REPORT THIS TO THE DEVELOPERS")
        raise ValueError

    def reset_iterator(self, chromosome: str, start: int):
        # changes chromosome and gets new iterator using .bai index file
        self.chromosome = chromosome
        self.reads_iterator: IteratorRowRegion = self.BAM_handle.fetch(f"{self.chromosome_prefix}{chromosome}",
                                                                       start=start, multiple_iterators=False)
        self.last_unmapped_read = self.get_next_mapped_read()
        self.last_extracted_reads = []

    def backtrack_reads(self, start: int, end: int) -> List[AlignedSegment]:
        # get reads that map from the last query
        new_last_extracted_reads = []
        ret = []
        for read in self.last_extracted_reads:
            if read.reference_start + 1 <= start and end <= read.reference_end:
                ret.append(read)
            elif end <= read.reference_end: # could theoretically map to the next locus if it's start and end are later
                new_last_extracted_reads.append(read)
        self.last_extracted_reads = new_last_extracted_reads
        return ret

    @staticmethod
    def simple_filter(read: AlignedSegment) -> bool:
        return not (read.flag & FLAG_OPTIONS.SECONDARY_ALG or read.flag & FLAG_OPTIONS.POOR_QUALITY
                    or read.flag & FLAG_OPTIONS.DUPLICATE_READ or read.flag & FLAG_OPTIONS.SUPPLEMENTARY_ALG)

    def get_next_mapped_read(self) -> AlignedSegment:
        cur_read = next(self.reads_iterator, None)
        while cur_read is not None:
            if cur_read.reference_end is None:  #  unaligned
                cur_read = next(self.reads_iterator, None)
            else:
                return cur_read
        return None

    def add_all_mapped(self, mapped_reads: List[AlignedSegment], cur_read: AlignedSegment, start: int, end: int) -> List[AlignedSegment]:
        # this function is called once we have found the first mapped read and we now add all other mapped reads
        while cur_read.reference_start + 1 <= start:
            if end <= cur_read.reference_end and self.simple_filter(cur_read):
                mapped_reads.append(cur_read)
            cur_read = self.get_next_mapped_read()
            if cur_read is None:  # reads iterator is exhausted
                break
        return self.remember_return(mapped_reads, cur_read)

    @staticmethod
    def strip_chromosome(chromosome: str):
        """ return chromosome as only its name(ex. chr16 -> 16, chrX -> X) """
        if chromosome.upper().find("X") != -1:
            return "X"
        elif chromosome.upper().find("Y") != -1:
            return "Y"
        elif chromosome.upper().find("M") != -1:
            return "M"
        else:
            numeric_filter = filter(str.isdigit, chromosome)
            return "".join(numeric_filter)

    def remember_return(self, mapped_reads: List[AlignedSegment], last_read: AlignedSegment) -> List[AlignedSegment]:
        #  logs reads for use in for next reads fetch and returns them
        self.last_unmapped_read = last_read
        self.last_extracted_reads+=mapped_reads
        return mapped_reads

    def get_reads(self, chromosome, start: int, end: int) -> List[AlignedSegment]:
        chromosome = self.strip_chromosome(chromosome)
        if chromosome != self.chromosome or self.last_unmapped_read is None or abs(start - self.last_unmapped_read.reference_start) > 7500:
            self.reset_iterator(chromosome, start)
        mapped_reads = self.backtrack_reads(start, end)
        cur_read = self.last_unmapped_read
        while cur_read is not None:
            if start < cur_read.reference_start + 1:
                return self.remember_return(mapped_reads, cur_read)
            elif cur_read.reference_start + 1 <= start and end <= cur_read.reference_end: # from pysam documentation: 'reference_end points to one past the last aligned residue'. Hence, no plus 1 even tho its 0 indexed and phobos is 1 indexed
                return self.add_all_mapped(mapped_reads, cur_read, start, end)
            else:
                cur_read = self.get_next_mapped_read()
        return self.remember_return(mapped_reads, cur_read)
