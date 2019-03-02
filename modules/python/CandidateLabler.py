"""
REUSED FROM HERE: https://github.com/google/deepvariant/blob/r0.7/deepvariant/labeler/haplotype_labeler.py
// REUSE DECLARATION: THIS SCRIPT HEAVILY REUSES THE METHODS USED IN DEEPVARIANT https://github.com/google/deepvariant/
// LICENSE INCLUDED IN: third_party/deepVariant.LICENSE
"""
from build import HELEN
from operator import itemgetter
import collections
import itertools
import heapq

DEBUG_IT = True
"""
            possible combos:
            gt1       gt2     Candidate validated?
            --------------------------------------
            hom       hom     no
            het       hom     yes
            het       het     yes
            hom_alt   hom     yes
            hom       None    no
            het       None    yes
            hom_alt   None    yes

            impossible combos:
            gt1       gt2     Candidate validated?
            --------------------------------------
            hom_alt   hom_alt NA
            het       hom_alt NA
            hom       hom_alt NA
            hom_alt   het     NA
            hom       het     NA
            None      None    NA

"""
MAX_GROUP_SIZE = 8
MAX_SEPARATION = 30

GENOTYPE_DICT = {"Hom": 0, "Het": 1, "Hom_alt": 2}
VCF_OFFSET = 1
PLOIDY = 2

# DEBUG_FREQUENCIES = False
DEBUG_PRINT_ALL = False

# Candidate data indexes
CHR_NAME = 0
START = 1
STOP = 2
REF_INDEX = 3
ALT1 = 4
ALT2 = 5
ALT1_TYPE = 6
ALT2_TYPE = 7

# Positional vcf indexes
SNP, IN, DEL = 0, 1, 2
SNP_CANDIDATE, IN_CANDIDATE, DEL_CANDIDATE = 1, 2, 3
SNP_DEL = 3

# VCF record indexes
REF, ALT, GT = 0, 1, 2

# Genotype codes
HOM, HET, HOM_ALT = 0, 1, 2

_CANDIDATE_MARKER = 'candidate'
_TRUTH_MARKER = 'truth'
_VariantToGroup = collections.namedtuple('_VariantToGroup',
                                         ['start', 'type', 'variant'])


def _variant_genotypes(variants, missing_genotypes_default=(-1, -1)):
    return [v[4] for v in variants]


def n_zeroes(l):
    return sum(1 for x in l if x == 0)


class HaplotypeMatch(object):
    def __init__(self, haplotypes, candidates, candidate_genotypes, truths, truth_genotypes):

        self.haplotypes = sorted(haplotypes)
        self.candidates = candidates
        self.truths = truths
        self.candidate_genotypes = candidate_genotypes
        self.truth_genotypes = truth_genotypes

        # Computed on-demand.
        self._n_false_positives = None
        self._n_false_negatives = None

    def __str__(self):
        return ('HaplotypeMatch(haplotypes={}, false_negatives={}, '
                'false_positives={} true_positives={} match_metrics={}, '
                'variant_gts={}, true_gts={})').format(
            self.haplotypes, self.n_false_negatives, self.n_false_positives,
            self.n_true_positives, self.match_metrics,
            self.candidate_genotypes, self.truth_genotypes)

    __repr__ = __str__

    @property
    def original_truth_genotypes(self):
        return _variant_genotypes(self.truths)

    @property
    def match_metrics(self):
        """Quality of this match. Lower scores are better.
        Returns:
          tuple[int] where all elements are >= 0: The tuple is suitable for sorting
          matches, so that sorted(matches, key=lambda x: x.match_metrics) will rank
          matches so that the best option is first.
        """
        return (self.n_false_negatives, self.n_false_positives,
                self.n_true_positives)

    @property
    def n_true_positives(self):
        """Gets the number of candidates whose matched genotype is not (0, 0).
        Since the candidates don't have expected genotypes, we can only count each
        site instead of each genotype. So this is the number of candidates whose
        matched genotype is not (0, 0).
        Returns:
          int >= 0.
        """
        return len(self.candidate_genotypes) - self.n_false_positives

    @property
    def n_false_positives(self):
        """Gets the number of candidates whose matched genotype is (0, 0).
        Since the candidates don't have expected genotypes, we can only count each
        site instead of each genotype. So this is the number of candidates whose
        matched genotype is (0, 0).
        Returns:
          int >= 0.
        """
        if self._n_false_positives is None:
            self._n_false_positives = sum(
                sum(gt) == 0 for gt in self.candidate_genotypes)
        return self._n_false_positives

    @property
    def n_false_negatives(self):
        """Gets the number of missed true genotypes.
        This is the sum of missed non-ref genotypes over all truth variants. So if
        we have a matched truth genotype of (0, 1) and the true genotype is (1, 1),
        then we have 1 FN. If the matched genotype were (0, 0), we'd have 2 FNs.
        Returns:
          int >= 0.
        """
        if self._n_false_negatives is None:
            self._n_false_negatives = sum(
                n_zeroes(assigned_gt) - n_zeroes(original_gt)
                for original_gt, assigned_gt in zip(self.original_truth_genotypes,
                                                    self.truth_genotypes))
        return self._n_false_negatives

    def candidates_with_assigned_genotypes(self):
        with_gts = []
        for variant, gt in zip(self.candidates, self.candidate_genotypes):
            v = (variant[0], variant[1], variant[2], variant[3], gt, variant[5], variant[6])
            with_gts.append(v)
        return with_gts


class RefCache:
    def __init__(self, ref_seq, ref_start, ref_end):
        self.ref_seq = ref_seq
        self.ref_start = ref_start
        self.ref_end = ref_end

    def get_seq(self, start_pos, end_pos):
        start_index = start_pos - self.ref_start
        end_index = end_pos - self.ref_start
        return self.ref_seq[start_index:end_index]


class ImpossibleHaplotype(Exception):
    """Indicates that an impossible haplotype configuration has been observed."""
    pass


class CandidateLabeler:
    def __init__(self, fasta_file_path, vcf_file_path):
        """
        Initialize candidateLabeler object
        """
        self.vcf_offset = -1
        self.delete_char = '*'
        self.fasta_file_path = fasta_file_path
        self.vcf_file_path = vcf_file_path

    def get_reference_sequence(self, candidate_group, truth_group):
        all_variants = candidate_group + truth_group
        contig = all_variants[0][0]
        start = min(v[1] for v in all_variants)
        end = max(v[2] for v in all_variants)
        start = start
        end = end
        fasta_handler = FRIDAY.FASTA_handler(self.fasta_file_path)
        ref_seq = fasta_handler.get_reference_sequence(contig, start, end)

        return RefCache(ref_seq, start, end)

    def all_possible_genotypes(self, gt):
        alts = set(gt) - {0}
        return {(0, 0), tuple(gt)} | {(0, alt) for alt in alts}

    def get_genotypes_for_truth(self, truth_set):
        gts = [v[4] for v in truth_set]
        return [self.all_possible_genotypes(gt) for gt in gts]

    def get_genotypes_for_candidate(self, num_alts):
        for j in range(num_alts + 1):
            for i in range(j + 1):
                yield (i, j)

    @staticmethod
    def is_overlapping(range_a, range_b):
        return min(range_a[1], range_b[1]) - max(range_a[0], range_b[0]) >= 0

    def split_independent_variants(self, variants_and_genotypes):
        overlaps = [variants_and_genotypes[0]]
        for i in range(1, len(variants_and_genotypes)):
            vgi = variants_and_genotypes[i][0]
            if any(self.is_overlapping((vg[0][1], vg[0][2]), (vgi[1], vgi[2])) for vg in overlaps):
                overlaps.append(variants_and_genotypes[i])
            else:
                return overlaps, variants_and_genotypes[i:]
        return overlaps, []

    def _allele_from_index(self, variant, allele_index):
        alleles = variant[3]
        return alleles[allele_index]

    def build_haplotype(self, variants, allele_indices, ref, ref_start, ref_end):
        parts = []
        position = ref_start
        for variant, allele_index in zip(variants, allele_indices):
            if variant[1] < position:
                if allele_index != 0:
                    return None
            else:
                ref_prefix = ref.get_seq(position, variant[1])
                # print("Build haplotype")
                # print("pos: ", position, variant[1], variant[:-1])
                # print("allele index", allele_index)
                # print("ref prefix: ", ref_prefix)
                allele = self._allele_from_index(variant, allele_index)
                if allele_index == 0:
                    allele = allele[0]
                    position = variant[1] + 1
                else:
                    position = variant[2]
                # print("allele: ", allele)
                parts.append(ref_prefix + allele)

        # We have some bases left to add between the position of our last variant
        # and the ref_end, so append those now.
        if position < ref_end:
            parts.append(ref.get_seq(position, ref_end))

        return ''.join(parts)

    def phased_genotypes_to_haplotypes(self, variants_and_genotypes, start, ref):
        genotypes_to_haplotypes = {}
        genotypes = [vg[1] for vg in variants_and_genotypes]
        variants = [vg[0] for vg in variants_and_genotypes]
        all_haploid_genotypes = sorted(set(itertools.product(*genotypes)))
        end = max(v[2] for v in variants)
        # print("------------")
        # print("Genotype to haplotype\n", variants_and_genotypes, start)
        for phased in all_haploid_genotypes:
            haplotype = self.build_haplotype(variants, phased, ref, start, end)
            if haplotype:
                genotypes_to_haplotypes[phased] = haplotype
        # print(genotypes_to_haplotypes)
        # print("------------")
        return genotypes_to_haplotypes, end

    def extend_haplotypes(self, prefix_haplotypes_list, haplotypes):
        for prefix_haplotypes in prefix_haplotypes_list:
            if len(prefix_haplotypes) == 1:
                f, = prefix_haplotypes
                yield {f + h for h in haplotypes}
            else:
                f1, f2 = prefix_haplotypes
                if len(haplotypes) == 1:
                    h, = haplotypes
                    yield {f1 + h, f2 + h}
                else:
                    h1, h2 = haplotypes
                    yield {f1 + h1, f2 + h2}
                    yield {f1 + h2, f2 + h1}

    def all_diploid_haplotypes(self, variants_and_genotypes, genotypes2haplotype):
        """Returns all diploid haplotypes for variants given their genotypes."""

        def complement_haploid_genotype(haploid_genotype, genotypes):
            assert len(haploid_genotype) == len(genotypes)
            return tuple(
                g1[1] if hg1 == g1[0] and len(g1) == 2 else g1[0]
                for hg1, g1 in zip(haploid_genotype, genotypes)
            )
        genotypes = [vg[1] for vg in variants_and_genotypes]
        generated_already = set()
        for haploid_genotype, haplotype in genotypes2haplotype.items():
            complement = complement_haploid_genotype(haploid_genotype, genotypes)
            complement_haplotype = genotypes2haplotype.get(complement, None)
            if complement_haplotype is not None and complement not in generated_already:
                generated_already.add(haploid_genotype)
                yield {haplotype, complement_haplotype}

    def create_haplotypes_recursive(self, variants_and_genotypes, ref, last_pos):
        if not variants_and_genotypes:
            yield {ref.get_seq(last_pos, ref.ref_end)} if last_pos != ref.ref_end else {''}
        else:
            group, remaining = self.split_independent_variants(variants_and_genotypes)
            group_haplotypes, next_pos = self.phased_genotypes_to_haplotypes(group, last_pos, ref)
            prefix_haplotypes = list(self.all_diploid_haplotypes(group, group_haplotypes))
            if not prefix_haplotypes:
                raise ImpossibleHaplotype
            for haplotypes in self.create_haplotypes_recursive(remaining, ref, next_pos):
                for result in self.extend_haplotypes(prefix_haplotypes, haplotypes):
                    yield result

    def create_haplotypes(self, variants_and_genotypes, ref, last_pos):
        try:
            for r in self.create_haplotypes_recursive(variants_and_genotypes, ref, last_pos):
                yield r
        except ImpossibleHaplotype:
            pass

    def create_truth_haplotype(self, truth_set, ref):
        all_possible_genotypes = self.get_genotypes_for_truth(truth_set)

        for genotypes in itertools.product(*all_possible_genotypes):
            combined = [(v, g) for v, g in zip(truth_set, genotypes)]
            for haplotypes in self.create_haplotypes(combined, ref, ref.ref_start):
                yield haplotypes, genotypes

    def get_genotypes_for_candidate_set(self, candidate_set):
        return [{gt for gt in self.get_genotypes_for_candidate(len(v[3]) - 1)} for v in candidate_set]

    def create_candidate_haplotype(self, candidate_set, ref):
        all_possible_genotypes = self.get_genotypes_for_candidate_set(candidate_set)
        # print("......................................................")
        # print("CANDIDATE GENOTYPES:", all_possible_genotypes)

        for genotypes in itertools.product(*all_possible_genotypes):
            combined = [(v, g) for v, g in zip(candidate_set, genotypes)]
            # print("COMBINED: ", combined)
            for haplotypes in self.create_haplotypes(combined, ref, ref.ref_start):
                yield haplotypes, genotypes

    def group_variants(self, candidate_records, chromosome_name, pos_start, pos_end):
        def _include_in_variant_group(group, group_variant):
            if not group:
                return True
            n_of_type = sum(1 for g in group if g.type == group_variant.type)
            if n_of_type >= MAX_GROUP_SIZE:
                return False
            else:
                return any(
                    group_variant.variant[1] - g.variant[2] + 1 <= MAX_SEPARATION
                    for g in group)

        def to_grouped_variants(variants, candidate_type):
            """Converts a Variant proto to a _VariantToGroup tuple."""
            return [_VariantToGroup(v[1], candidate_type, v) for v in variants]

        vcf_handler = FRIDAY.VCF_handler(self.vcf_file_path)
        truth_records = vcf_handler.get_vcf_records(chromosome_name, pos_start - 1, pos_end + 1)
        converted_records = []
        for t_rec in truth_records:
            single_rec = (t_rec.chromosome_name, t_rec.start_pos, t_rec.end_pos,
                          tuple(t_rec.alleles), tuple(t_rec.genotype), True)
            converted_records.append(single_rec)
        truth_records = converted_records

        # print("CANDIDATES:")
        # for candidate in candidate_records:
        #     print(candidate)
        # print("--------------------------")
        # print("TRUE")
        # for true_record in truth_records:
        #     print(true_record)
        # print("--------------------------")

        groupable_variants = heapq.merge(
            to_grouped_variants(candidate_records, _CANDIDATE_MARKER),
            to_grouped_variants(truth_records, _TRUTH_MARKER))

        # split the candidate group
        groups = []
        current_group = []
        for record in groupable_variants:
            can_add = _include_in_variant_group(current_group, record)

            if can_add:
                current_group.append(record)
            else:
                groups.append(current_group)

                current_group = [record]

        if current_group:
            groups.append(current_group)

        separated_groups = []
        for group in groups:
            candidates = [record.variant for record in group if record.type == _CANDIDATE_MARKER]
            truths = [record.variant for record in group if record.type == _TRUTH_MARKER]
            separated_groups.append((candidates, truths))

        return separated_groups

    @staticmethod
    def solve_multiple_alts(alts, ref):
        type1, type2 = alts[0][1], alts[1][1]
        alt1, alt2 = alts[0][0], alts[1][0]
        if type1 == DEL_CANDIDATE and type2 == DEL_CANDIDATE:
            if len(alt2) > len(alt1):
                return alt2, alt2[0] + alt2[len(alt1):], ref
            else:
                return alt1, ref, alt1[0] + alt1[len(alt2):]
        elif type1 == IN_CANDIDATE and type2 == IN_CANDIDATE:
            return ref, alt1, alt2
        elif type1 == DEL_CANDIDATE or type2 == DEL_CANDIDATE:
            if type1 == DEL_CANDIDATE and type2 == IN_CANDIDATE:
                return alt1, ref, alt2 + alt1[1:]
            elif type1 == IN_CANDIDATE and type2 == DEL_CANDIDATE:
                return alt2, alt1 + alt2[1:], ref
            elif type1 == DEL_CANDIDATE and type2 == SNP_CANDIDATE:
                return alt1, ref, alt2 + alt1[1:]
            elif type1 == SNP_CANDIDATE and type2 == DEL_CANDIDATE:
                return alt2, alt1 + alt2[1:], ref
            elif type1 == DEL_CANDIDATE:
                return alt1, ref, alt2
            elif type2 == DEL_CANDIDATE:
                return alt2, alt1, ref
        else:
            return ref, alt1, alt2

    @staticmethod
    def solve_single_alt(alts, ref):
        # print(alts)
        alt1, alt_type = alts
        if alt_type == DEL_CANDIDATE:
            return alt1, ref
        return ref, alt1

    def get_proper_candidate(self, candidate):
        if candidate.alt2_type == 0:
            allele_tuple = self.solve_single_alt((candidate.alt1, candidate.alt1_type), candidate.ref)
        else:
            allele_tuple = self.solve_multiple_alts(((candidate.alt1, candidate.alt1_type),
                                                     (candidate.alt2, candidate.alt2_type)), candidate.ref)

        return candidate.chromosome_name, candidate.pos, candidate.pos + len(allele_tuple[0]), \
               allele_tuple, (0, 0), False, candidate

    def duplicate_haplotypes(self, haplotypes_and_genotypes):
        haplotypes_and_genotypes = list(haplotypes_and_genotypes)
        return [
            (vh1, vg1)
            for i, (vh1, vg1) in enumerate(haplotypes_and_genotypes)
            if not any(vh1 == vh2 for (vh2, _) in haplotypes_and_genotypes[i + 1:])
        ]

    def select_best_haplotype_match(self, all_matches):
        sorted_matches = sorted(all_matches, key=lambda x: x.match_metrics)
        best = sorted_matches[0]
        equivalents = [
            f for f in all_matches if f.match_metrics == best.match_metrics
        ]

        return equivalents[0]

    def find_best_haplotype(self, candidate_group, truth_group, ref):
        truth_haplotypes = self.duplicate_haplotypes(self.create_truth_haplotype(truth_group, ref))
        candidate_haplotypes = list(self.create_candidate_haplotype(candidate_group, ref))
        # if DEBUG_IT:
        #     print("TRUTH HAPLOTYPES:")
        #     for hap in truth_haplotypes:
        #         print(hap)
        #     print("CANDIDATE HAPLOTYPES:")
        #     for hap in candidate_haplotypes:
        #         print(hap)

        found = []
        for vh, vgt in candidate_haplotypes:
            for th, tgt in truth_haplotypes:
                if th == vh:
                    found.append(
                        HaplotypeMatch(
                            haplotypes=th,
                            candidates=candidate_group,
                            candidate_genotypes=vgt,
                            truths=truth_group,
                            truth_genotypes=tgt))
        del candidate_haplotypes, truth_haplotypes
        if not found:
            return None
        else:
            return self.select_best_haplotype_match(found)

    def get_labeled_candidates(self, chromosome_name, pos_start, pos_end, candidate_sites):
        """
        Label candidates given variants from a VCF
        :param positional_vcf: IN/DEL/SNPs separated into VCF records, expanded into a 1-to-1 ref_pos:variant allele
        :param candidate_sites: Candidates
        :return: List of labeled candidate sites
        """
        # list of all labeled candidates
        candidate_records = []

        # for candidate in candidate_sites:
        #     candidate.print()

        for candidate in candidate_sites:
            candidate_record = self.get_proper_candidate(candidate)
            candidate_records.append(candidate_record)

        del candidate_sites

        # sort candidates based on positions
        candidate_records = sorted(candidate_records, key=itemgetter(1))
        groups = self.group_variants(candidate_records, chromosome_name, pos_start, pos_end)

        all_labeled_candidates = []
        for candidate_group, truth_group in groups:
            # if DEBUG_IT:
            #     print("-------------START-----------")
            #     print("Candidate group:")
            #     for candidate in candidate_group:
            #         print(candidate[:-1])
            #     print("########################")
            #     print("Truth group:")
            #     for candidate in truth_group:
            #         print(candidate[:-1])
            #     print("........................")
            if not candidate_group:
                continue
            if not truth_group:
                for labeled_candidate in candidate_group:
                    labeled_candidate[6].set_genotype(labeled_candidate[4])
                    # candidate_with_gts = labeled_candidate[6].set_genotype(labeled_candidate[4])
                    all_labeled_candidates.append(labeled_candidate[6])
                continue

            ref = self.get_reference_sequence(candidate_group, truth_group)
            # if DEBUG_IT:
            #     print("REF:")
            #     print(ref.ref_seq, ref.ref_start, ref.ref_end)
            #     print(".......................")
            labeled_set = self.find_best_haplotype(candidate_group, truth_group, ref)

            if labeled_set is None:
                raise ValueError('Failed to assign labels for variants', ref.ref_start, ref.ref_end)

            # if DEBUG_IT:
            #     print("LABELED CANDIDATES:")
            for labeled_candidate in labeled_set.candidates_with_assigned_genotypes():
                labeled_candidate[6].set_genotype(labeled_candidate[4])
                # if DEBUG_IT:
                #     print(candidate_with_gts)
                all_labeled_candidates.append(labeled_candidate[6])

            # if DEBUG_IT:
            #     print("------------------------")
            del labeled_set, ref

        # all_labeled_candidates = sorted(all_labeled_candidates, key=itemgetter(1))
        # for labeled_candidate in all_labeled_candidates:
        #     labeled_candidate.print()
        # labeled_candidate.print()
        # print(all_labeled_candidates)
        return all_labeled_candidates
