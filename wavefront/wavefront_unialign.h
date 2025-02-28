/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef WAVEFRONT_UNIALIGN_H_
#define WAVEFRONT_UNIALIGN_H_

#include "wavefront_aligner.h"

/*
 * Resize
 */
void wavefront_unialign_resize(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const bool reverse_sequences);

/*
 * Initialize alignment
 */
void wavefront_unialign_initialize_wavefronts(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length);
void wavefront_unialign_init(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end);

/*
 * Classic WF-Alignment (Unidirectional)
 */
int wavefront_unialign(
    wavefront_aligner_t* const wf_aligner);

/*
 * Display
 */
void wavefront_unialign_print_status(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner,
    const int current_score);
void wavefront_align_unidirectional_cleanup(
    wavefront_aligner_t* const wf_aligner);
wf_offset_t accumulate_converge_idx(wf_offset_t a,wf_offset_t b);
wf_offset_t test_converged_score(wavefront_aligner_t* const wf_aligner, int score);
void wavefront_unialign_terminate(
    wavefront_aligner_t* const wf_aligner,
    const int score);
#endif /* WAVEFRONT_UNIALIGN_H_ */

