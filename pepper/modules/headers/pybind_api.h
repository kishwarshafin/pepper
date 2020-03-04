//
// Created by Kishwar Shafin on 10/18/18.
//

#ifndef PEPPER_PYBIND_API_H
#define PEPPER_PYBIND_API_H

#include "dataio/fasta_handler.h"
#include "dataio/bam_handler.h"
#include "realignment/simple_aligner.h"
#include "pileup_summary/summary_generator.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

PYBIND11_MODULE(PEPPER, m) {
        py::class_<SummaryGenerator>(m, "SummaryGenerator")
            .def(py::init<const string &, const string &, long long &, long long &>())
            .def_readwrite("genomic_pos", &SummaryGenerator::genomic_pos)
            .def_readwrite("labels", &SummaryGenerator::labels)
            .def_readwrite("image", &SummaryGenerator::image)
            .def_readwrite("bad_label_positions", &SummaryGenerator::bad_label_positions)
            .def("generate_train_summary", &SummaryGenerator::generate_train_summary)
            .def("generate_summary", &SummaryGenerator::generate_summary);

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
            .def("align_reads_to_reference", &ReadAligner::align_reads_to_reference);

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
}
#endif //PEPPER_PYBIND_API_H
