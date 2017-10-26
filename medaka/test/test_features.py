import os
import numpy as np
import unittest
from medaka import features

from medaka.inference import load_feature_file, save_feature_file, run_prediction
from medaka.inference import load_encoding, generate_samples, load_model_hdf
from medaka.tview import load_pileup, rechunk, generate_pileup_chunks
from collections import Counter, OrderedDict


__bam_fp__ = os.path.join(os.path.dirname(__file__), 'data', 'test_reads.bam')
__model_fp__ = os.path.join(os.path.dirname(__file__), 'data', 'test_model.hdf5')
__encoding_fp__ = os.path.join(os.path.dirname(__file__), 'data', 'test_encodings.json')
__ref_fasta__ = os.path.join(os.path.dirname(__file__), 'data', 'draft_ref.fasta')
__ref_name__ = 'Consensus_Consensus_Consensus_Consensus_utg000001l'


def gen_overlaps(data, field='features'):

    for s1, s2 in zip(data[0:-1], data[1:]):

        start_ind1 = np.searchsorted(s1.positions, s2.positions[0])
        end_ind2 = np.searchsorted(s2.positions, s1.positions[-1], side='right')

        pos1_ovl = s1.positions[start_ind1:]
        pos2_ovl = s2.positions[0:end_ind2]

        field1_ovl = getattr(s1, field)[start_ind1:]
        field2_ovl = getattr(s2, field)[0:end_ind2]

        yield pos1_ovl, pos2_ovl, field1_ovl, field2_ovl


class TviewTest(unittest.TestCase):

    def setUp(self):

        self.tview_overlap = 1000
        self.rechunk_overlap = 200

        self.tview_gen = generate_pileup_chunks((__bam_fp__,), __ref_fasta__, __ref_name__, start=50000, end=250000, 
                                                overlap=self.tview_overlap, chunk_len=25000)


    def test_overlapping_counts_same_within_tview_chunk(self):
        """After rechunking a single tview chunk, check counts features in overlapping regions are identical"""
        tview_chunk = next(self.tview_gen)
        small_chunks = rechunk((tview_chunk,), chunk_size=1000, overlap=self.rechunk_overlap)
        data = generate_samples(small_chunks, coverage_filter=iter,
                                      feature_func=features.counts)
        save_feature_file('rechunk_features.npy', data)
        for pos1, pos2, feat1, feat2 in gen_overlaps(data):
            assert len(pos1) == len(pos2)
            assert np.all(pos1 == pos2)
            assert len(pos1) >= self.rechunk_overlap
            assert np.all(feat1 == feat2)


    def test_overlapping_counts_in_adjacent_tview_chunks(self):
        """Check counts features in overlapping regions are identical"""

        feature_file = 'tview_chunk_features.npy'
        reuse_features = False
        if reuse_features and os.path.exists(feature_file):
            data = load_feature_file(feature_file)
        else:
            data = generate_samples(self.tview_gen, coverage_filter=iter,
                                          feature_func=features.counts)
            save_feature_file(feature_file, data)
        counts = Counter()
        for pos1, pos2, feat1, feat2 in gen_overlaps(data):
            assert len(pos1) == len(pos2)
            assert len(np.unique(pos1['major'])) == self.tview_overlap
            assert len(np.unique(pos2['major'])) == self.tview_overlap
            assert pos2['major'][-1] - pos2['major'][0] == self.tview_overlap - 1
            assert pos1['major'][-1] - pos1['major'][0] == self.tview_overlap - 1
            # these two could fail if some reads are included in one tview chunk 
            # but not the other but that does not seem to happen (at least for our test set)
            assert np.all(pos1 == pos2)
            counts.update([np.max(np.abs(feat1 - feat2))])
        # it seems this is not true at two tview chunk borders, where we have a diff of 1 base
        # due to an short segment of a read at the start of the second chunk not being printed by tview
        # check 0 is most common value
        assert counts[0] == max(counts.values())
        # check we have max count diff of 1, would ideally be zero
        assert max(counts.keys()) == 1 


    def test_inference(self):
        """Check predictions in overlapping regions"""

        big_chunk = next(self.tview_gen)
        small_chunks = rechunk((big_chunk,), chunk_size=1000, overlap=self.rechunk_overlap)

        data = generate_samples(small_chunks, coverage_filter=iter,
                                feature_func=features.counts)

        for pos1, pos2, feat1, feat2 in gen_overlaps(data):
            assert len(pos1) == len(pos2)
            assert np.all(pos1 == pos2)
            assert len(pos1) >= self.rechunk_overlap
            assert np.all(feat1 == feat2)

        save_feature_file('inference_features.hdf', data)
        pred_file = 'inference_predictions.hdf'
        model, encoding = load_model_hdf(__model_fp__, encoding_json=__encoding_fp__)
        pred_data = run_prediction(data, model, encoding, 
                       output_file='rechunk_basecalls.fasta', batch_size=1500,
                       predictions_file=pred_file)

        # it seems that in some cases there may be differences when the most probable and second
        # most probable have almost equal probability 
        counts = Counter()
        differences = []
        for pos1, pos2, pred1, pred2 in gen_overlaps(pred_data, field='label_probs'):
            assert len(pos1) == len(pos2)
            assert np.all(pos1 == pos2)
            assert len(pos1) >= self.rechunk_overlap
            best1 = np.argmax(pred1, -1)
            best2 = np.argmax(pred2, -1)
            same = np.all(best1 == best2)
            counts.update([same])
            if not same:
                differ_inds = np.where(best1 != best2)[0]
                for ind in differ_inds:
                    row = pred1[ind]
                    sort_inds = np.argsort(row)
                    best = sort_inds[-1]
                    second = sort_inds[-2]
                    coverage = np.sum(feat1[ind])
                    differences.append((pos1[ind], ind, best, row[best], second, row[second], coverage, feat1[ind]))
        # we should have only 1 difference in 1 of the 49 small chunks
        if counts[False] > 1:
            raise AssertionError("Expected <= 1 of the small chunks to differ in some way, got {}".format(counts))
        if len(differences) > 1:
            msg = "position, index, best, p(best), second_best, p(second_best), coverage, counts_feature\n"
            for i in differences:
                msg += str(i) + '\n'
            raise AssertionError("More differences found than expected:\n {}".format(msg))
        assert counts[True] == 49  # check something worked


if __name__ == '__main__':
    unittest.main()
