import pysam

"""
This class handles bam files using pysam API.
"""


class BamHandler:
    """
    Handles bam files using pysam API
    """
    def __init__(self, bam_file_path):
        """
        create AlignmentFile object given file path to a bam file
        :param bam_file_path: full path to a bam file
        """
        try:
            self.bamFile = pysam.AlignmentFile(bam_file_path, "rb")
        except:
            raise IOError("BAM FILE READ ERROR")

    def get_pileupcolumns_aligned_to_a_region(self, chromosome_name, start, end):
        """
        Return a AlignmentFile.pileup object given a site
        :param contig: Contig [ex. chr3]
        :param pos: Position [ex 100001]
        :return: pysam.AlignmentFile.pileup object
        """
        # get pileup columns
        pileup_columns = self.bamFile.pileup(chromosome_name, start, end, truncate=True)
        # return pileup columns
        return pileup_columns

    def get_pileupcolumns_aligned_to_a_site(self, contig, pos):
        """
        Return a AlignmentFile.pileup object given a site
        :param contig: Contig [ex. chr3]
        :param pos: Position [ex 100001]
        :return: pysam.AlignmentFile.pileup object
        """
        # get pileup columns
        pileup_columns = self.bamFile.pileup(contig, pos, pos+1)
        # return pileup columns
        return pileup_columns

    def get_reads(self, chromosome_name, start, stop):
        """
        Return reads that map to a given site
        :param chromosome_name: Chromosome name. Ex: chr3
        :param start: Site start in the chromosome
        :param stop: Site end in the chromosome
        :return: Reads that align to that site
        """
        return self.bamFile.fetch(chromosome_name, start, stop)

    def get_header_sq(self):
        return self.bamFile.header['SQ']
