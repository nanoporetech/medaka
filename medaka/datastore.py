from collections import Counter, defaultdict, OrderedDict
import h5py
import numpy as np
import logging
import yaml

from medaka.common import Sample, decoding


class DataStore(object):
    """Class to read/write to a data file"""
    _sample_path_ = 'samples'
    _groups_ = ('medaka_features_kwargs', 'medaka_model_kwargs',
                'medaka_label_decoding', 'medaka_feature_decoding',
                'medaka_label_counts')

    def __init__(self, filename, mode='r'):

        self.filename = filename
        self.mode = mode

        self.samples_written = set()
        self.fh = None

        self.logger = logging.getLogger(__package__)
        self.logger.name = 'DataStore'

        self._meta = None


    def __enter__(self):

        self.fh = h5py.File(self.filename, self.mode)

        return self


    def __exit__(self, *args):
        if self.mode != 'r' and self._meta is not None:
            self._write_metadata(self.meta)
        self.fh.close()


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

        :param sample: `Sample` object.
        """
        # If we are writing samples, count labels and store them in meta
        if 'medaka_label_counts' not in self.meta:
            self.meta['medaka_label_counts'] = Counter()

        if sample.name not in self.samples_written:
            self.samples_written.add(sample.name)
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
        else:
            self.logger.debug('Not writing {} as it is present already'.format(sample.name))


    def load_sample(self, key):
        """Load `Sample` object from HDF5

        :param key: str, sample name.
        :returns: `Sample` object.
        """
        s = {}
        for field in Sample._fields:
            pth = '{}/{}/{}'.format(self._sample_path_, key, field)
            if pth in self.fh:
                s[field] = self.fh[pth][()]
                if isinstance(s[field], np.ndarray) and isinstance(s[field][0], type(b'')):
                    s[field] = np.char.decode(s[field])
            else:
                s[field] = None
        return Sample(**s)


    def log_counts(self):
        """Log label counts"""

        h = lambda l: (decoding[l[0]], l[1]) if type(l) == tuple else l
        self.logger.info("Label counts:\n{}".format('\n'.join(
            ['{}: {}'.format(h(label), count) for label, count in self.meta['medaka_label_counts'].items()]
        )))


    def _write_metadata(self, data):
        """Save a data structure to file within a yml str."""
        for group, d in data.items():
            if group in self.fh:
                del self.fh[group]
            self.fh[group] = yaml.dump(d)


    def _load_metadata(self, groups=None):
        """Load meta data"""
        if groups is None:
            groups = self._groups_
        return {g: yaml.load(self.fh[g][()]) for g in groups if g in self.fh}

    @property
    def sample_keys(self):
        """Return list of sample keys"""
        return tuple(self.fh[self._sample_path_].keys()) if self._sample_path_ in self.fh else ()

    @property
    def n_samples(self):
        """Return number of samples"""
        return len(self.sample_keys)


class DataIndex(object):
    """Class to index and serve samples from one or more `DataFiles`"""

    def __init__(self, filenames):

        self.filenames = filenames

        with DataStore(filenames[0]) as ds:
            self.meta = ds.meta

        c_grp = 'medaka_label_counts'
        if c_grp in self.meta:
            self.meta[c_grp] = Counter()

        self.samples = []
        for f in filenames:
            with DataStore(f) as ds:
                self.samples.extend([(s, f) for s in ds.sample_keys])
                if c_grp in self.meta:
                    self.meta[c_grp].update(ds._load_metadata(groups=(c_grp,))[c_grp])

        self._index = None

        self.logger = logging.getLogger(__package__)
        self.logger.name = 'DataIndex'


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
            d = Sample.decode_sample_name(key)
            if d is not None:
                d['key'] = key
                d['filename'] = f
                ref_names[d['ref_name']].append(d)

        # sort dicts so that refs are in order and within a ref, chunks are in order
        ref_names_ordered = OrderedDict()
        for ref_name in sorted(ref_names.keys()):
            sorter = lambda x: float(x['start'])
            ref_names[ref_name].sort(key=sorter)
            ref_names_ordered[ref_name] = ref_names[ref_name]

        return ref_names_ordered


    def yield_from_feature_files(self, ref_names=None, samples=None):
        """Yield `Sample` objects from one or more feature files.

        :ref_names: iterable of str, only process these references.
        :samples: iterable of sample names to yield (in order in which they are supplied).
        :yields: `Sample` objects.
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
