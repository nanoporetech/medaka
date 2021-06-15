#ifndef _MEDAKA_RNNVARS_H
#define _MEDAKA_RNNVARS_H

/** Calculate which columns of an RNN output matrix should be considered
 * variant.
 *
 * @param minor pileup minor indices (non-zero elements count insertion).
 * @param reference reference sequence including gaps for insertion pileup
 *     columns.
 * @param prediction sequence predicted by RNN including gaps.
 * @param bool (out) whether pileup column is considered part of a variant run.
 *     It should be zero-initialised (no variants)
 * @param len the length of all input and output arrays.  @returns void
 *
 * The function produces the following: - A reference pileup column is variant
 * if it contains a substitution or deletion.  - An insertion pileup column is
 * variant if any pileup column (reference or insertion) associated with the
 * same reference position is variant.
 *
 * Note that a reference position is not variant if simply the training
 * insertion columns are variant. i.e. the function is not attempting to create
 * boundaries in the vein of a VCF record -- only annotating columns containing
 * variant bases.
 *
 */
void variant_columns(size_t* minor, wchar_t* reference, wchar_t* prediction, bool* out, size_t len);

#endif
