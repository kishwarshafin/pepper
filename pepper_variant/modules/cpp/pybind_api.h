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
#include "candidate_finder_hp.h"
#include "region_summary.h"
#include "region_summary_hp.h"
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
            .def_readwrite("chunked_type_labels", &RegionalImageSummary::chunked_type_labels)
            .def_readwrite("chunked_ids", &RegionalImageSummary::chunked_ids);

    // Generate regional summary
       py::class_<RegionalSummaryGenerator>(m, "RegionalSummaryGenerator")
            .def(py::init<string &, long long &, long long &, string &>())
            .def_readwrite("max_observed_insert", &RegionalSummaryGenerator::max_observed_insert)
            .def_readwrite("cumulative_observed_insert", &RegionalSummaryGenerator::cumulative_observed_insert)
            .def_readwrite("total_observered_insert_bases", &RegionalSummaryGenerator::total_observered_insert_bases)
            .def("generate_summary", &RegionalSummaryGenerator::generate_summary)
            .def("generate_labels", &RegionalSummaryGenerator::generate_labels)
            .def("generate_max_insert_summary", &RegionalSummaryGenerator::generate_max_insert_summary);

        py::class_<RegionalSummaryGeneratorHP>(m, "RegionalSummaryGeneratorHP")
                .def(py::init<string &, long long &, long long &, string &>())
                .def_readwrite("max_observed_insert", &RegionalSummaryGeneratorHP::max_observed_insert)
                .def_readwrite("cumulative_observed_insert", &RegionalSummaryGeneratorHP::cumulative_observed_insert)
                .def_readwrite("total_observered_insert_bases", &RegionalSummaryGeneratorHP::total_observered_insert_bases)
                .def("generate_summary", &RegionalSummaryGeneratorHP::generate_summary)
                .def("generate_labels", &RegionalSummaryGeneratorHP::generate_labels)
                .def("generate_max_insert_summary", &RegionalSummaryGeneratorHP::generate_max_insert_summary);

        py::class_<CandidateImageSummary>(m, "CandidateImageSummary")
            .def(py::init<>())
            .def(py::init<string &, long long &, int &, vector<string> &, vector<int> &, vector<vector<int> > &, int &, int & >())
            .def_readwrite("contig", &CandidateImageSummary::contig)
            .def_readwrite("position", &CandidateImageSummary::position)
            .def_readwrite("depth", &CandidateImageSummary::depth)
            .def_readwrite("candidates", &CandidateImageSummary::candidates)
            .def_readwrite("candidate_frequency", &CandidateImageSummary::candidate_frequency)
            .def_readwrite("image_matrix", &CandidateImageSummary::image_matrix)
            .def_readwrite("base_label", &CandidateImageSummary::base_label)
            .def_readwrite("type_label", &CandidateImageSummary::type_label)
            .def(py::pickle(
                    [](const CandidateImageSummary &p) { // __getstate__
                        /* Return a tuple that fully encodes the state of the object */
                        return py::make_tuple(p.contig, p.position, p.depth, p.candidates, p.candidate_frequency, p.image_matrix, p.base_label, p.type_label);
                        },
                        [](py::tuple t) { // __setstate__
                            if (t.size() != 8)
                                throw std::runtime_error("Invalid state!");

                            /* Create a new C++ instance */
                            CandidateImageSummary p(t[0].cast<string>(), t[1].cast<long long>(), t[2].cast<int>(), t[3].cast< vector<string> >(), t[4].cast< vector<int> >(), t[5].cast<vector<vector<int> > >(), t[6].cast<int>(), t[7].cast<int>() );

                            /* Assign any additional state */
                            //dp.setExtra(t[1].cast<int>());

                            return p;
                        }
                ));

        // HP image generator
        py::class_<CandidateImageSummaryHP>(m, "CandidateImageSummaryHP")
                .def(py::init<>())
                .def(py::init<string &, long long &, int &, vector<string> &, vector<int> &, vector<vector<int> > &, int &, int & >())
                .def_readwrite("contig", &CandidateImageSummaryHP::contig)
                .def_readwrite("position", &CandidateImageSummaryHP::position)
                .def_readwrite("depth", &CandidateImageSummaryHP::depth)
                .def_readwrite("candidates", &CandidateImageSummaryHP::candidates)
                .def_readwrite("candidate_frequency", &CandidateImageSummaryHP::candidate_frequency)
                .def_readwrite("image_matrix", &CandidateImageSummaryHP::image_matrix)
                .def_readwrite("base_label", &CandidateImageSummaryHP::base_label)
                .def_readwrite("type_label", &CandidateImageSummaryHP::type_label)
                .def(py::pickle(
                        [](const CandidateImageSummaryHP &p) { // __getstate__
                            /* Return a tuple that fully encodes the state of the object */
                            return py::make_tuple(p.contig, p.position, p.depth, p.candidates, p.candidate_frequency, p.image_matrix, p.base_label, p.type_label);
                        },
                        [](py::tuple t) { // __setstate__
                            if (t.size() != 8)
                                throw std::runtime_error("Invalid state!");

                            /* Create a new C++ instance */
                            CandidateImageSummaryHP p(t[0].cast<string>(), t[1].cast<long long>(), t[2].cast<int>(), t[3].cast< vector<string> >(), t[4].cast< vector<int> >(), t[5].cast<vector<vector<int> > >(), t[6].cast<int>(), t[7].cast<int>() );

                            /* Assign any additional state */
                            //dp.setExtra(t[1].cast<int>());

                            return p;
                        }
                ));
        py::class_<CandidateImagePrediction>(m, "CandidateImagePrediction")
            .def(py::init<>())
            .def(py::init<string &, long long &, int &, vector<string> &, vector<int> &, vector<float> &, vector<float> &>())
            .def_readwrite("contig", &CandidateImagePrediction::contig)
            .def_readwrite("position", &CandidateImagePrediction::position)
            .def_readwrite("depth", &CandidateImagePrediction::depth)
            .def_readwrite("candidates", &CandidateImagePrediction::candidates)
            .def_readwrite("candidate_frequency", &CandidateImagePrediction::candidate_frequency)
            .def_readwrite("prediction_base", &CandidateImagePrediction::prediction_base)
            .def_readwrite("prediction_type", &CandidateImagePrediction::prediction_type)
            .def(py::pickle(
                    [](const CandidateImagePrediction &p) { // __getstate__
                        /* Return a tuple that fully encodes the state of the object */
                        return py::make_tuple(p.contig, p.position, p.depth, p.candidates, p.candidate_frequency, p.prediction_base, p.prediction_type);
                    },
                    [](py::tuple t) { // __setstate__
                        if (t.size() != 7)
                            throw std::runtime_error("Invalid state CandidateImagePrediction!");

                        /* Create a new C++ instance */
                        CandidateImagePrediction p(t[0].cast<string>(), t[1].cast<long long>(), t[2].cast<int>(), t[3].cast<vector<string> >(), t[4].cast<vector<int> >(), t[5].cast<vector<float> >(), t[6].cast< vector<float> >());

                        /* Assign any additional state */
                        //dp.setExtra(t[1].cast<int>());

                        return p;
                    }
            ));


        // data structure for methylation record
        py::class_<type_truth_record>(m, "type_truth_record")
                .def(py::init<string &, long long &, long long &, string &, string &>())
                .def_readwrite("contig", &type_truth_record::contig)
                .def_readwrite("pos_start", &type_truth_record::pos_start)
                .def_readwrite("pos_end", &type_truth_record::pos_end)
                .def_readwrite("ref", &type_truth_record::ref)
                .def_readwrite("alt", &type_truth_record::alt);

        // data structure for methylation record
        py::class_<type_truth_recordHP>(m, "type_truth_recordHP")
                .def(py::init<string &, long long &, long long &, string &, string &>())
                .def_readwrite("contig", &type_truth_recordHP::contig)
                .def_readwrite("pos_start", &type_truth_recordHP::pos_start)
                .def_readwrite("pos_end", &type_truth_recordHP::pos_end)
                .def_readwrite("ref", &type_truth_recordHP::ref)
                .def_readwrite("alt", &type_truth_recordHP::alt);

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
            .def(py::init<const string &, const string &, long long &, long long&, long long &, long long &>())
            .def("find_candidates_consensus", &CandidateFinder::find_candidates_consensus)
            .def("find_candidates", &CandidateFinder::find_candidates);

        py::class_<CandidateFinderHP>(m, "CandidateFinderHP")
            .def(py::init<const string &, const string &, long long &, long long&, long long&, long long&>())
            .def("find_candidates", &CandidateFinderHP::find_candidates);

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
            .def_readwrite("allele_probability", &Candidate::allele_probability)
            .def_readwrite("genotype_probability", &Candidate::genotype_probability)
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
