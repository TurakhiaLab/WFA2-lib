extern "C"{
#include "utils/commons.h"
#include "wavefront/wavefront_align.h"
}
#include <fstream>

int main(int argc, char** argv){
      wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine;
  attributes.affine_penalties.match = 0;
  attributes.affine_penalties.mismatch = 1;
  attributes.affine_penalties.gap_opening = 1;
  attributes.affine_penalties.gap_extension = 1;
  attributes.use_tile=true;
    // Set heuristic wf-adaptive
  attributes.heuristic.strategy = wf_heuristic_wfadaptive;
  attributes.heuristic.min_wavefront_length = 0;
  attributes.heuristic.max_distance_threshold = 10;
  attributes.heuristic.steps_between_cutoffs = 1;
  //Semi-global
      attributes.alignment_form.span = alignment_endsfree;
      attributes.alignment_form.pattern_begin_free = 0;
    attributes.alignment_form.pattern_end_free = 1;
    attributes.alignment_form.text_begin_free = 0;
    attributes.alignment_form.text_end_free = 1;
  attributes.memory_mode = wavefront_memory_med;

  // Initialize Wavefront Aligner
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  wf_aligner->marking_score=16;
        
    int aligned=0;
    int failed=0;
    std::fstream file (argv[1], std::fstream::in);
    while (true) {
        std::string ref,seq_name,query;
        std::getline(file, seq_name);
        std::getline(file, query);
        std::getline(file, seq_name);
        if(seq_name==""){
            break;
        }
        std::getline(file, ref);
        char* cigar;
        int cigar_len;
        Effi_Stats_t stats;
        wavefront_tile(wf_aligner,&cigar,&cigar_len,
        query.data(),
        query.size(),
        ref.data(),
        ref.size(),
        &stats);
  cigar_t cigar_inst;
  cigar_inst.begin_offset=0;
  cigar_inst.end_offset=cigar_len-1;
  cigar_inst.operations=cigar;
  fprintf(stderr,"  PATTERN  %s\n",query.c_str());
  fprintf(stderr,"  TEXT     %s\n",ref.c_str());
  fprintf(stderr,"  SCORE (RE)COMPUTED %d\n",
      cigar_score_gap_affine(&cigar_inst,&attributes.affine_penalties));
  cigar_print_pretty(stderr,
      query.c_str(),query.size(),ref.c_str(),ref.size(),
      &cigar_inst,wf_aligner->mm_allocator);
    free(cigar);
    }
}