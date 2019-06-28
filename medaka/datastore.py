from collections import Counter, defaultdict, OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import yaml

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py

import medaka.common


class DataStore(object):
    """Class to read/write to a data file"""
    _sample_path_ = 'samples'
    _groups_ = ('medaka_features_kwargs', 'medaka_model_kwargs', 'medaka_model_name',
                'medaka_label_decoding', 'medaka_feature_decoding',
                'medaka_label_counts', 'medaka_samples')

    def __init__(self, filename, mode='r', verify_on_close=True):

        self.filename = filename
        self.mode = mode
        self.verify_on_close = verify_on_close

        self._sample_keys = set()
        self.fh = None

        self.logger = medaka.common.get_named_logger('DataStore')

        self._meta = None


    def __enter__(self):

        self.fh = h5py.File(self.filename, self.mode)

        return self


    def __exit__(self, *args):

        if self.mode != 'r':
            if self.verify_on_close:
                self._verify_()  # verify data before saving meta
            else:
                self.logger.debug("Skipping validation on close.")
            self._write_metadata(self.meta)
        self.fh.close()


    def _verify_(self):
        self.logger.debug("Verifying data.")
        self.fh.flush()
        fh = h5py.File(self.filename, 'r')
        # find the union of all present fields and remove and samples from the
        # index which don't contain all of these fields
        all_fields = set()

        for key in self.sample_keys:
            self.logger.debug("First round verify {}.".format(key))
            # if key is not in the file, remove it from the index
            grp = '{}/{}'.format(self._sample_path_, key)
            if grp not in fh:
                self.meta['medaka_samples'].remove(key)
                self.logger.debug("Removing sample {} as grp {} is not present.".format(key, grp))
                continue
            all_fields.update(fh[grp].keys())

        for key in self.sample_keys:
            self.logger.debug("Second round verify {}.".format(key))
            for field in all_fields:
                path = '{}/{}/{}'.format(self._sample_path_, key, field)
                if path not in fh:
                    self.meta['medaka_samples'].remove(key)
                    self.logger.debug("Removing sample {} as field {} is not present.".format(key, path))
                    break


    @property
    def meta(self):
        if self._meta is None:
            self._meta = self._load_metadata()
        return self._meta


    def update_meta(self, meta):
        """Update metadata"""
        self._meta = self.meta
        self._meta.update(meta)


    def write_sample(self, sample):
        """Write sample to hdf, ensuring a sample is not written twice and maintaining
        a count of labels seen.

        :param sample: `medaka.common.Sample` object.
        """
        # count labels and store them in meta
        if 'medaka_label_counts' not in self.meta:
            self.meta['medaka_label_counts'] = Counter()
        # Store sample index in meta
        if 'medaka_samples' not in self.meta:
            self.meta['medaka_samples'] = set()

        if not any([isinstance(getattr(sample, field), np.ndarray) for field in sample._fields]):
            self.logger.debug('Not writing sample as it has no data.')
        elif sample.name not in self.meta['medaka_samples']:
            for field in sample._fields:
                if getattr(sample, field) is not None:
                    data = getattr(sample, field)
                    if isinstance(data, np.ndarray) and isinstance(data[0], np.unicode):
                        data = np.char.encode(data)
                    self.fh['{}/{}/{}'.format(self._sample_path_, sample.name, field)] = data

            if sample.labels is not None:
                if len(sample.labels.dtype) == 2:  # RLE-encoded
                    self.meta['medaka_label_counts'].update([tuple(l) for l in sample.labels])
                else:
                    self.meta['medaka_label_counts'].update(sample.labels)
            # Do this last so we only add this sample to the index if we have
            # gotten this far
            self.meta['medaka_samples'].add(sample.name)
        else:
            self.logger.debug('Not writing {} as it is present already'.format(sample.name))


    def load_sample(self, key):
        """Load `medaka.common.Sample` object from HDF5

        :param key: str, sample name.
        :returns: `medaka.common.Sample` object.
        """
        s = {}
        for field in medaka.common.Sample._fields:
            pth = '{}/{}/{}'.format(self._sample_path_, key, field)
            if pth in self.fh:
                s[field] = self.fh[pth][()]
                if isinstance(s[field], np.ndarray) and isinstance(s[field][0], type(b'')):
                    s[field] = np.char.decode(s[field])
            else:
                s[field] = None
        return medaka.common.Sample(**s)


    def log_counts(self):
        """Log label counts"""

        h = lambda l: (medaka.common.decoding[l[0]], l[1]) if type(l) == tuple else l
        self.logger.info("Label counts:\n{}".format('\n'.join(
            ['{}: {}'.format(h(label), count) for label, count in self.meta['medaka_label_counts'].items()]
        )))


    def _write_metadata(self, data):
        """Save a data structure to file within a yml str."""
        self.logger.debug("Writing metadata.")
        for group, d in data.items():
            if group in self.fh:
                del self.fh[group]
            self.fh[group] = yaml.dump(d)


    def _load_metadata(self, groups=None):
        """Load meta data"""
        if groups is None:
            groups = self._groups_
        return {g: yaml.unsafe_load(self.fh[g][()]) for g in groups if g in self.fh}


    @property
    def sample_keys(self):
        """Return tuple of sample keys"""

        return tuple(self.meta['medaka_samples']) if 'medaka_samples' in self.meta else tuple()


    def _find_samples(self):
        """Find samples in a file and update meta with a list of samples."""
        self.update_meta({'medaka_samples': set(self.fh[self._sample_path_].keys()) if self._sample_path_ in self.fh else {}})


    @property
    def n_samples(self):
        """Return number of samples"""
        return len(self.sample_keys)


class DataIndex(object):
    """Class to index and serve samples from one or more `DataFiles`"""

    def __init__(self, filenames, threads=4):

        self.logger = medaka.common.get_named_logger('DataIndex')

        self.filenames = filenames

        with DataStore(filenames[0]) as ds:
            self.logger.debug('Loading meta from {}'.format(filenames[0]))
            self.meta = ds.meta

        c_grp = 'medaka_label_counts'
        if c_grp in self.meta:
            self.meta[c_grp] = Counter()

        if 'medaka_samples' in self.meta:
            del self.meta['medaka_samples']

        self.samples = []

        with ProcessPoolExecutor(threads) as executor:
            future_to_f = {executor.submit(DataIndex._load_meta, f): f for f in filenames}
            for i, future in enumerate(as_completed(future_to_f), 1):
                f = future_to_f[future]
                try:
                    meta = future.result()
                    if 'medaka_samples' in meta:
                        self.samples.extend([(s, f) for s in meta['medaka_samples']])
                    else:
                        self.logger.info('Could not find samples in {}'.format(f))

                    self.meta[c_grp].update(meta[c_grp])
                except Exception as exc:
                    self.logger.info('Could not load meta from {}'.format(f))
                else:
                    self.logger.info('Loaded sample-index from {}/{} ({:.2%}) of feature files.'.format(i, len(filenames), i / len(filenames)))

        # make order of samples independent of order in which tasks complete
        self.samples.sort()

        self._index = None

    @staticmethod
    def _load_meta(f):
        with DataStore(f) as ds:
            meta = ds.meta
            #medaka.common.get_named_logger('Load_meta').debug('Done {}'.format(f))
            return meta


    @property
    def index(self):
        self._index = self._get_sorted_index() if self._index is None else self._index
        return self._index


    def _get_sorted_index(self):
        """Get index of samples indexed by reference and ordered by start pos.

        :returns: {ref_name: [sample dicts sorted by start]}
        """

        ref_names = defaultdict(list)

        for key, f in self.samples:
            d = medaka.common.Sample.decode_sample_name(key)
            if d is not None:
                d['key'] = key
                d['filename'] = f
                ref_names[d['ref_name']].append(d)

        # sort dicts so that refs are in order and within a ref, chunks are in order
        ref_names_ordered = OrderedDict()

        get_major_minor = lambda x: tuple((int(i) for i in x.split('.')))
        # sort by start and -end so that if we have two samples with the same
        # start but differrent end points, the longest sample comes first
        sorter = lambda x: (get_major_minor(x['start']) + tuple((-i for i in get_major_minor(x['end']))))
        for ref_name in sorted(ref_names.keys()):
            ref_names[ref_name].sort(key=sorter)
            ref_names_ordered[ref_name] = ref_names[ref_name]

        return ref_names_ordered


    def yield_from_feature_files(self, ref_names=None, samples=None):
        """Yield `medaka.common.Sample` objects from one or more feature files.

        :ref_names: iterable of str, only process these references.
        :samples: iterable of sample names to yield (in order in which they are supplied).
        :yields: `medaka.common.Sample` objects.
        """

        if samples is not None:
            # yield samples in the order they are asked for
            for sample, fname in samples:
                yield DataStore(fname).load_sample(sample)
        else:
            # yield samples sorted by ref_name and start
            if ref_names is None:
                ref_names = sorted(self.index.keys())
            for ref_name in ref_names:
                for d in self.index[ref_name]:
                    with DataStore(d['filename']) as ds:
                        yield ds.load_sample(d['key'])
