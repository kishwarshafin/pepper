//
// Created by Kishwar Shafin on 10/18/18.
//

#ifndef PEPPER_VARIANT_PYBIND_API_H
#define PEPPER_VARIANT_PYBIND_API_H

#include "cigar.h"
#include "read.h"
#include "sequence.h"
#include "fasta_handler.h"
#include "bam_handler.h"
#include "summary_generator.h"
#include "candidate_finder.h"
#include "region_summary.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

PYBIND11_MODULE(PEPPER_VARIANT, m) {
        py::class_<ImageSummary>(m, "ImageSummary")
            .def_readwrite("images", &ImageSummary::images)
            .def_readwrite("positions", &ImageSummary::positions)
            .def_readwrite("refs", &ImageSummary::refs)
            .def_readwrite("labels", &ImageSummary::labels)
            .def_readwrite("chunk_ids", &ImageSummary::chunk_ids);

        // Generate summary
        py::class_<SummaryGenerator>(m, "SummaryGenerator")
            .def(py::init<const string &, const string &, long long &, long long &>())
            .def_readwrite("genomic_pos", &SummaryGenerator::genomic_pos)
            .def_readwrite("ref_image", &SummaryGenerator::ref_image)
            .def_readwrite("labels", &SummaryGenerator::labels)
            .def_readwrite("image", &SummaryGenerator::image)
            .def_readwrite("bad_label_positions", &SummaryGenerator::bad_label_positions)
            .def_readwrite("longest_insert_count", &SummaryGenerator::longest_insert_count)
            .def("generate_train_summary", &SummaryGenerator::generate_train_summary)
            .def("chunk_image", &SummaryGenerator::chunk_image)
            .def("chunk_image_train", &SummaryGenerator::chunk_image_train)
            .def("generate_summary", &SummaryGenerator::generate_summary);

        //
        py::class_<RegionalImageSummary>(m, "RegionalImageSummary")
            .def_readwrite("chunked_image_matrix", &RegionalImageSummary::chunked_image_matrix)
            .def_readwrite("chunked_positions", &RegionalImageSummary::chunked_positions)
            .def_readwrite("chunked_index", &RegionalImageSummary::chunked_index)
            .def_readwrite("chunked_labels", &RegionalImageSummary::chunked_labels)
            .def_readwrite("chunked_ids", &RegionalImageSummary::chunked_ids);

       // Generate regional summary
       py::class_<RegionalSummaryGenerator>(m, "RegionalSummaryGenerator")
            .def(py::init<long long &, long long &, string &>())
            .def_readwrite("max_observed_insert", &RegionalSummaryGenerator::max_observed_insert)
            .def_readwrite("cumulative_observed_insert", &RegionalSummaryGenerator::cumulative_observed_insert)
            .def_readwrite("total_observered_insert_bases", &RegionalSummaryGenerator::total_observered_insert_bases)
            .def("generate_summary", &RegionalSummaryGenerator::generate_summary)
            .def("generate_labels", &RegionalSummaryGenerator::generate_labels)
            .def("generate_max_insert_summary", &RegionalSummaryGenerator::generate_max_insert_summary);

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

        // Candidate finder
        py::class_<CandidateFinder>(m, "CandidateFinder")
            .def(py::init<const string &, const string &, long long &, long long&, long long&, long long&>())
            .def("find_candidates", &CandidateFinder::find_candidates);

        py::class_<CandidateAllele>(m, "CandidateAllele")
            .def(py::init<>())
            .def_readwrite("ref", &CandidateAllele::ref)
            .def_readwrite("alt", &CandidateAllele::alt)
            .def_readwrite("alt_type", &CandidateAllele::alt_type);

        py::class_<Candidate>(m, "Candidate")
            .def(py::init<>())
            .def("print", &Candidate::print)
            .def("set_genotype", &Candidate::set_genotype)
            .def_readwrite("pos_start", &Candidate::pos)
            .def_readwrite("pos_end", &Candidate::pos_end)
            .def_readwrite("read_support", &Candidate::read_support)
            .def_readwrite("alt_prob", &Candidate::alt_prob)
            .def_readwrite("alt_prob_h1", &Candidate::alt_prob_h1)
            .def_readwrite("alt_prob_h2", &Candidate::alt_prob_h2)
            .def_readwrite("non_ref_prob", &Candidate::non_ref_prob)
            .def_readwrite("depth", &Candidate::depth)
            .def_readwrite("genotype", &Candidate::genotype)
            .def_readwrite("allele", &Candidate::allele);

        py::class_<PositionalCandidateRecord>(m, "PositionalCandidateRecord")
            .def(py::init<>())
            .def(py::init<const string &, long long &, long long&, const int &>())
            .def(py::self < py::self)
            .def_readwrite("chromosome_name", &PositionalCandidateRecord::chromosome_name)
            .def_readwrite("pos_start", &PositionalCandidateRecord::pos_start)
            .def_readwrite("pos_end", &PositionalCandidateRecord::pos_end)
            .def_readwrite("candidates", &PositionalCandidateRecord::candidates);
}
#endif //PEPPER_VARIANT_PYBIND_API_H
