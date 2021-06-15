
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

/** Calculate which columns of an RNN output matrix should be considered
 * variant.
 *
 * @param minor pileup minor indices (non-zero elements count insertion).
 * @param reference reference sequence including gaps for insertion pileup
 * columns.
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
void variant_columns(size_t* minor, wchar_t* reference, wchar_t* prediction, bool* out, size_t len) {
    bool is_var = reference[0] != prediction[0];  // assume start on major
    size_t insert_length = 0;
    out[0] = is_var;
    for (size_t i = 1; i < len; ++i) {
        if (minor[i] == 0) {
            // start of new reference position
            if (is_var) {
                // if we saw any vars in an insert run, set all inserts to true
                for (size_t j = i - insert_length; j < i; ++j) {
                    out[j] = true;
                }
            }
            is_var = reference[i] != prediction[i];
            out[i] = is_var;
            insert_length = 0;
        } else {
            insert_length += 1;
            is_var = is_var || reference[i] != prediction[i];
        }
    }
    // set any remaining inserts
    if (is_var) {
        for (size_t j = len - insert_length; j <= len - 1; ++j) {
            out[j] = true;
        }
    }
}
