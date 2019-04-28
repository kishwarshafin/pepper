from build import HELEN
label_decoder = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: ''}


def valid_base(base):
    if base == 'A' or base == 'a' or \
            base == 'C' or base == 'c' or \
            base == 'G' or base == 'g' or \
            base == 't' or base == 'T':
        return True
    return False


class Haplotype2Variant:
    def __init__(self, fasta_file_path, contig, region_start, region_end):
        self.contig = contig
        self.region_start = region_start
        self.region_end = region_end
        self.fasta_handler = HELEN.FASTA_handler(fasta_file_path)

    def generate_records_from_haplotypes(self, predictions):
        ref_start = self.region_start
        ref_end = self.region_end
        ref_sequence = self.fasta_handler.get_reference_sequence(self.contig,
                                                                 ref_start,
                                                                 ref_end)

        variant_set = list()
        for pos in range(ref_start, ref_end - 1):
            indx = 0
            ref_base = ref_sequence[pos-ref_start]
            haplotype_1_base = label_decoder[predictions[(pos, indx, 1)]]
            haplotype_2_base = label_decoder[predictions[(pos, indx, 2)]]

            if not valid_base(ref_base):
                continue

            alternate_alleles = []
            genotype = [0, 0]
            if haplotype_1_base != ref_base and valid_base(haplotype_1_base):
                alternate_alleles.append(haplotype_1_base)
            if haplotype_2_base != ref_base and valid_base(haplotype_2_base):
                alternate_alleles.append(haplotype_2_base)

            if len(alternate_alleles) == 1:
                genotype = [0, 1]
            elif len(alternate_alleles) == 2:
                if alternate_alleles[0] == alternate_alleles[1]:
                    alternate_alleles = [alternate_alleles[0]]
                    genotype = [1, 1]
                else:
                    genotype = [1, 2]

            if len(alternate_alleles) > 0:
                variant_set.append((self.contig, pos, pos+1, ref_base, alternate_alleles, genotype))

        return variant_set
