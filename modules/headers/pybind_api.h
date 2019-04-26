//
// Created by Kishwar Shafin on 10/18/18.
//

#ifndef HELEN_PYBIND_API_H
#define HELEN_PYBIND_API_H

#include "dataio/fasta_handler.h"
#include "dataio/bam_handler.h"
#include "dataio/vcf_handler.h"
#include "local_reassembly/active_region_finder.h"
#include "local_reassembly/debruijn_graph.h"
#include "local_reassembly/aligner.h"
#include "candidate_finding/candidate_finder.h"
#include "image_generator/image_generator.h"
#include "pileup_summary/summary_generator.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

PYBIND11_MODULE(HELEN, m) {
        py::class_<PileupImage>(m, "PileupImage")
            .def(py::init<>())
            .def_readwrite("chromosome_name", &PileupImage::chromosome_name)
            .def_readwrite("start_pos", &PileupImage::start_pos)
            .def_readwrite("end_pos", &PileupImage::end_pos)
            .def_readwrite("image", &PileupImage::image)
            .def_readwrite("label", &PileupImage::label);

        py::class_<ImageGenerator>(m, "ImageGenerator")
            .def(py::init<const string &, const string &, long long &, long long&,
        map<long long, PositionalCandidateRecord> & >())
            .def("create_window_pileups", &ImageGenerator::create_window_pileups)
            .def("set_positional_vcf", &ImageGenerator::set_positional_vcf);

        py::class_<SummaryGenerator>(m, "SummaryGenerator")
            .def(py::init<const string &, const string &, long long &, long long &>())
            .def_readwrite("genomic_pos", &SummaryGenerator::genomic_pos)
            .def_readwrite("labels", &SummaryGenerator::labels)
            .def_readwrite("image", &SummaryGenerator::image)
            .def_readwrite("bad_label_positions", &SummaryGenerator::bad_label_positions)
            .def("generate_train_summary", &SummaryGenerator::generate_train_summary)
            .def("generate_summary", &SummaryGenerator::generate_summary);

        py::class_<PositionalCandidateRecord>(m, "PositionalCandidateRecord")
            .def(py::init<>())
            .def(py::init<const string &, long long &, long long&, const string &,
                 const string &, const string &, const int &, const int &>())
            .def("print", &PositionalCandidateRecord::print)
            .def("set_genotype", &PositionalCandidateRecord::set_genotype)
            .def("get_candidate_record", &PositionalCandidateRecord::get_candidate_record)
            .def_readwrite("chromosome_name", &PositionalCandidateRecord::chromosome_name)
            .def_readwrite("pos", &PositionalCandidateRecord::pos)
            .def_readwrite("pos_end", &PositionalCandidateRecord::pos_end)
            .def_readwrite("ref", &PositionalCandidateRecord::ref)
            .def_readwrite("alt1", &PositionalCandidateRecord::alt1)
            .def_readwrite("alt2", &PositionalCandidateRecord::alt2)
            .def_readwrite("alt1_type", &PositionalCandidateRecord::alt1_type)
            .def_readwrite("alt2_type", &PositionalCandidateRecord::alt2_type)
            .def(py::pickle(
                    [](const PositionalCandidateRecord &p) { // __getstate__
                        /* Return a tuple that fully encodes the state of the object */
                        return py::make_tuple(p.chromosome_name, p.pos, p.pos_end, p.ref, p.alt1, p.alt2, p.alt1_type, p.alt2_type);
                    },
                    [](py::tuple t) { // __setstate__
                        if (t.size() != 8)
                            throw std::runtime_error("Invalid state!");

                        /* Create a new C++ instance */
                        PositionalCandidateRecord p(t[0].cast<string>(), t[1].cast<long long>(), t[2].cast<long long>(), t[3].cast<string>(),
                                                    t[4].cast<string>(), t[5].cast<string>(), t[6].cast<int>(), t[7].cast<int>());

                        /* Assign any additional state */
                        //dp.setExtra(t[1].cast<int>());

                        return p;
                    }
            ));


        // Candidate finder
        py::class_<CandidateFinder>(m, "CandidateFinder")
            .def(py::init<const string &, const string &, long long &, long long&, long long&, long long&>())
            .def("find_candidates", &CandidateFinder::find_candidates);

        // Alignment CLASS
        py::class_<StripedSmithWaterman::Alignment>(m, "Alignment")
            .def(py::init<>())
            .def_readwrite("best_score", &StripedSmithWaterman::Alignment::sw_score)
            .def_readwrite("best_score2", &StripedSmithWaterman::Alignment::sw_score_next_best)
            .def_readwrite("reference_begin", &StripedSmithWaterman::Alignment::ref_begin)
            .def_readwrite("reference_end", &StripedSmithWaterman::Alignment::ref_end)
            .def_readwrite("query_begin", &StripedSmithWaterman::Alignment::query_begin)
            .def_readwrite("query_end", &StripedSmithWaterman::Alignment::query_end)
            .def_readwrite("ref_end_next_best", &StripedSmithWaterman::Alignment::ref_end_next_best)
            .def_readwrite("mismatches", &StripedSmithWaterman::Alignment::mismatches)
            .def_readwrite("cigar_string", &StripedSmithWaterman::Alignment::cigar_string)
            .def_readwrite("cigar", &StripedSmithWaterman::Alignment::cigar)
            .def("Clear", &StripedSmithWaterman::Alignment::Clear);

        // Filter Class
        py::class_<StripedSmithWaterman::Filter>(m, "Filter")
            .def_readwrite("report_begin_position", &StripedSmithWaterman::Filter::report_begin_position)
            .def_readwrite("report_cigar", &StripedSmithWaterman::Filter::report_cigar)
            .def_readwrite("score_filter", &StripedSmithWaterman::Filter::score_filter)
            .def_readwrite("distance_filter", &StripedSmithWaterman::Filter::distance_filter)
            .def(py::init<>())
            .def(py::init<const bool&, const bool&, const uint16_t&, const uint16_t&>());

        // Aligner Class
        py::class_<StripedSmithWaterman::Aligner>(m, "Aligner")
            .def(py::init<>())
            .def(py::init<const uint8_t&, const uint8_t&, const uint8_t&, const uint8_t&>())
            .def("SetReferenceSequence", &StripedSmithWaterman::Aligner::SetReferenceSequence)
            .def("Align_cpp", &StripedSmithWaterman::Aligner::Align_cpp);

        // SSW alignment class
        py::class_<ReadAligner>(m, "ReadAligner")
            .def(py::init<const int &, const int &, const string &>())
            .def("align_reads_to_reference", &ReadAligner::align_reads_to_reference)
            .def("align_reads", &ReadAligner::align_reads);

        // Debruijn graph class
        py::class_<DeBruijnGraph>(m, "DeBruijnGraph")
            .def(py::init<const long long &, const long long &>())
            .def_readwrite("current_hash_value", &DeBruijnGraph::current_hash_value)
            .def_readwrite("node_hash_int_to_str", &DeBruijnGraph::node_hash_int_to_str)
            .def_readwrite("good_nodes", &DeBruijnGraph::good_nodes)
            .def_readwrite("out_nodes", &DeBruijnGraph::out_nodes)
            .def_readwrite("edges", &DeBruijnGraph::edges)

            .def("generate_haplotypes", &DeBruijnGraph::generate_haplotypes)
            .def("find_min_k_from_ref", &DeBruijnGraph::find_min_k_from_ref);

        // data structure for sequence name and their length
        py::class_<ActiveRegionFinder>(m, "ActiveRegionFinder")
            .def(py::init<const string &, const string &, long long &, long long&>())
            .def("find_active_region", &ActiveRegionFinder::find_active_region);

        // data structure for sequence name and their length
        py::class_<type_sequence>(m, "type_sequence")
            .def_readwrite("sequence_length", &type_sequence::sequence_length)
            .def_readwrite("sequence_name", &type_sequence::sequence_name);

        // data structure for CIGAR operation
        py::class_<CigarOp>(m, "CigarOp")
            .def_readwrite("cigar_op", &CigarOp::operation)
            .def_readwrite("cigar_len", &CigarOp::length);

        // data structure for read attributes aka read flags
        py::class_<type_read_flags>(m, "type_read_flags")
            .def_readwrite("is_paired", &type_read_flags::is_paired)
            .def_readwrite("is_proper_pair", &type_read_flags::is_proper_pair)
            .def_readwrite("is_unmapped", &type_read_flags::is_unmapped)
            .def_readwrite("is_mate_unmapped", &type_read_flags::is_mate_unmapped)
            .def_readwrite("is_reverse", &type_read_flags::is_reverse)
            .def_readwrite("is_mate_is_reverse", &type_read_flags::is_mate_is_reverse)
            .def_readwrite("is_read1", &type_read_flags::is_read1)
            .def_readwrite("is_read2", &type_read_flags::is_read2)
            .def_readwrite("is_secondary", &type_read_flags::is_secondary)
            .def_readwrite("is_qc_failed", &type_read_flags::is_qc_failed)
            .def_readwrite("is_duplicate", &type_read_flags::is_duplicate)
            .def_readwrite("is_supplementary", &type_read_flags::is_supplementary)
            .def(py::init());

        // data structure for read
        py::class_<type_read>(m, "type_read")
            .def("set_read_id", &type_read::set_read_id)
            .def("__lt__", &type_read::operator<, py::is_operator())
            .def_readwrite("pos", &type_read::pos)
            .def_readwrite("pos_end", &type_read::pos_end)
            .def_readwrite("query_name", &type_read::query_name)
            .def_readwrite("read_id", &type_read::read_id)
            .def_readwrite("flags", &type_read::flags)
            .def_readwrite("hp_tag", &type_read::hp_tag)
            .def_readwrite("sequence", &type_read::sequence)
            .def_readwrite("cigar_tuples", &type_read::cigar_tuples)
            .def_readwrite("mapping_quality", &type_read::mapping_quality)
            .def_readwrite("base_qualities", &type_read::base_qualities)
            .def_readwrite("bad_indicies", &type_read::bad_indicies);

        // bam handler API
        py::class_<BAM_handler>(m, "BAM_handler")
            .def(py::init<const string &>())
            .def("get_chromosome_sequence_names", &BAM_handler::get_chromosome_sequence_names)
            .def("get_chromosome_sequence_names_with_length", &BAM_handler::get_chromosome_sequence_names_with_length)
            .def("get_reads", &BAM_handler::get_reads);

        // FASTA handler API
        py::class_<FASTA_handler>(m, "FASTA_handler")
            .def(py::init<const string &>())
            .def("get_reference_sequence", &FASTA_handler::get_reference_sequence)
            .def("get_chromosome_sequence_length", &FASTA_handler::get_chromosome_sequence_length)
            .def("get_chromosome_names", &FASTA_handler::get_chromosome_names);

        // VCF handler API
        py::class_<VCF_handler>(m, "VCF_handler")
            .def(py::init<const string &>())
            .def("get_vcf_records", &VCF_handler::get_vcf_records)
            .def("get_positional_vcf_records", &VCF_handler::get_positional_vcf_records);

        // VCF handler API
        py::class_<type_vcf_record>(m, "type_vcf_record")
            .def_readwrite("chromosome_name", &type_vcf_record::chromosome_name)
            .def_readwrite("start_pos", &type_vcf_record::start_pos)
            .def_readwrite("end_pos", &type_vcf_record::end_pos)
            .def_readwrite("id", &type_vcf_record::id)
            .def_readwrite("qual", &type_vcf_record::qual)
            .def_readwrite("is_filter_pass", &type_vcf_record::is_filter_pass)
            .def_readwrite("sample_name", &type_vcf_record::sample_name)
            .def_readwrite("genotype", &type_vcf_record::genotype)
            .def_readwrite("filters", &type_vcf_record::filters)
            .def_readwrite("alleles", &type_vcf_record::alleles);


        py::class_<type_alt_allele>(m, "type_alt_allele")
            .def_readonly("ref", &type_alt_allele::ref)
            .def_readonly("alt_allele", &type_alt_allele::alt_allele)
            .def_readonly("alt_type", &type_alt_allele::alt_type)
            .def("get_ref", &type_alt_allele::get_ref)
            .def("get_alt_allele", &type_alt_allele::get_alt_allele)
            .def("get_alt_type", &type_alt_allele::get_alt_type);

        py::class_<type_positional_vcf_record>(m, "type_positional_vcf_record")
            .def_readonly("chromosome_name", &type_positional_vcf_record::chromosome_name)
            .def_readonly("start_pos", &type_positional_vcf_record::start_pos)
            .def_readonly("end_pos", &type_positional_vcf_record::end_pos)
            .def_readonly("id", &type_positional_vcf_record::id)
            .def_readonly("qual", &type_positional_vcf_record::qual)
            .def_readonly("is_filter_pass", &type_positional_vcf_record::is_filter_pass)
            .def_readonly("sample_name", &type_positional_vcf_record::sample_name)
            .def_readonly("genotype", &type_positional_vcf_record::genotype)
            .def_readonly("filters", &type_positional_vcf_record::filters)
            .def_readonly("alt_allele", &type_positional_vcf_record::alt_allele);

}
#endif //HELEN_PYBIND_API_H
