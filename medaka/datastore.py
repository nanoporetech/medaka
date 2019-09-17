"""Storing of training and inference data to file."""
from collections import defaultdict, OrderedDict
from concurrent.futures import \
    as_completed, ProcessPoolExecutor, ThreadPoolExecutor
import pickle
import warnings

import numpy as np

import medaka.common

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class DataStore(object):
    """Read and write data to .hdf files."""

    _meta_group_ = 'meta'  # top level group for meta items
    _sample_path_ = 'samples/data'  # data group contains sample Datasets
    _sample_registry_path_ = 'samples/registry'  # set of sample keys

    def __init__(self, filename, mode='r'):
        """Initialize a datastore.

        :param filename: file to open.
        :param mode: file opening mode ('r', 'w', 'a').
        """
        self.filename = filename
        self.mode = mode

        self.logger = medaka.common.get_named_logger('DataStore')

        self.write_executor = ThreadPoolExecutor(1)
        self.write_futures = []

        self._sample_registry = None
        self.fh = h5py.File(self.filename, self.mode)

    def __enter__(self):
        """Create context for handling a datastore file."""
        return self

    def __exit__(self, *args):
        """Shutdown sample writer and close."""
        if self.mode != 'r':
            self.write_executor.shutdown(wait=True)
            self._write_sample_registry()
        self.close()

    def close(self):
        """Close filehandle of back-end file."""
        self.fh.close()

    def get_meta(self, key):
        """Load (deserialise) a meta data item.

        :param key: name of item to load.
        """
        path = '{}/{}'.format(self._meta_group_, key)
        try:
            return pickle.loads(self.fh[path][()])
        except Exception as e:
            self.logger.debug("Could not load {} from {}. {}.".format(
                key, self.filename, e))

    def set_meta(self, obj, key):
        """Store (serialize) a meta data item to file.

        :param obj: the object to serialise.
        :param key: the name of the object.
        """
        path = '{}/{}'.format(self._meta_group_, key)
        if path in self.fh:
            del self.fh[path]
        pickled_obj = np.string_(pickle.dumps(obj))
        self.fh[path] = pickled_obj
        self.fh.flush()

    def copy_meta(self, other):
        """Copy metadata to another file.

        :param other: filename of another file.
        """
        with DataStore(other, 'a') as other_ds:
            self.fh.copy(
                self._meta_group_, other_ds.fh,
                name=self._meta_group_)

    @property
    def sample_registry(self):
        """Return a set of samples stored within file."""
        self._initialise_sample_registry()  # is idempotent
        return self._sample_registry

    @property
    def n_samples(self):
        """Return the number of samples stored in file."""
        return len(self.sample_registry)

    def write_sample(self, sample):
        """Write sample to hdf.

        Checks are performed to ensure a sample is not written twice and
        a count of unique training labels seen is maintained.

        :param sample: `medaka.common.Sample` object.
        """
        contains_numpy_array = any(
            isinstance(getattr(sample, field), np.ndarray)
            for field in sample._fields)
        if not contains_numpy_array:
            self.logger.debug('Not writing sample as it has no data.')

        # if the sample does not already exist (according to sample registry)
        elif sample.name not in self.sample_registry:
            for field in sample._fields:
                # do not write None
                if getattr(sample, field) is not None:
                    data = getattr(sample, field)
                    # handle numpy array of unicode chars
                    if isinstance(data, np.ndarray) and \
                            isinstance(data[0], np.unicode):
                        data = np.char.encode(data)
                    location = '{}/{}/{}'.format(
                        self._sample_path_, sample.name, field)
                    self.write_futures.append(
                        self.write_executor.submit(
                            self._write_dataset, location, data))
            self._sample_registry.add(sample.name)
        else:
            self.logger.debug(
                'Not writing {} as present already'.format(
                    sample.name))

    def load_sample(self, key):
        """Load `medaka.common.Sample` object from file.

        :param key: str, sample name.
        :returns: `medaka.common.Sample` object.
        """
        s = dict()
        for field in medaka.common.Sample._fields:
            pth = '{}/{}/{}'.format(self._sample_path_, key, field)
            try:
                s[field] = self.fh[pth][()]
            except KeyError:
                s[field] = None
            else:
                # handle loading of bytestrings
                if isinstance(s[field], np.ndarray) and \
                        isinstance(s[field][0], type(b'')):
                    s[field] = np.char.decode(s[field])
        return medaka.common.Sample(**s)

    def _write_dataset(self, location, data):
        """Write data, compressing numpy arrays."""
        if isinstance(data, np.ndarray):
            self.fh.create_dataset(location, data=data,
                                   compression='gzip',
                                   compression_opts=1)
        else:
            self.fh[location] = data

    def _write_pickled(self, obj, path):
        """Write a pickled object to file."""
        if path in self.fh:
            del self.fh[path]
        pickled_obj = np.string_(pickle.dumps(obj))
        self.fh[path] = pickled_obj

    def _initialise_sample_registry(self):
        """Load sample registry from file if present else create."""
        if self._sample_registry is None:
            try:
                self._sample_registry = pickle.loads(
                    self.fh[self._sample_registry_path_][()])
            except KeyError:
                self._sample_registry = set()

    def _write_sample_registry(self):
        """Write sample registry."""
        self.logger.debug("Writing sample registry.")
        if self._sample_registry_path_ in self.fh:
            del self.fh[self._sample_registry_path_]
        self._write_pickled(self.sample_registry,
                            self._sample_registry_path_)


class DataIndex(object):
    """Index and serve samples from multiple `DataStore` compatible files."""

    def __init__(self, filenames, threads=4):
        """Intialize an index across a set of files.

        :param filenames: list of files to index.
        :param threads: number of threads to use for indexing.
        """
        self.logger = medaka.common.get_named_logger('DataIndex')
        self.filenames = filenames
        self.threads = threads
        self.n_files = len(self.filenames)
        self._index = None
        self._extract_sample_registries()
        self.metadata = self._load_metadata()

    def _extract_sample_registries(self):
        """."""
        self.samples = []

        with ProcessPoolExecutor(self.threads) as executor:
            future_to_fn = {
                executor.submit(DataIndex._load_sample_registry, fn): fn
                for fn in self.filenames}
            for i, future in enumerate(as_completed(future_to_fn), 1):
                fn = future_to_fn[future]
                try:
                    sample_registry = future.result()
                    self.samples.extend(
                        [(s, fn) for s in sample_registry])
                except Exception:
                    self.logger.info('No sample_registry in {}'.format(fn))
                else:
                    self.logger.info(
                        'Loaded {}/{} ({:.2f}%) sample files.'.format(
                            i, self.n_files,
                            i / self.n_files * 100))
        # make order of samples independent of order in which tasks complete
        self.samples.sort()

    @staticmethod
    def _load_sample_registry(f):
        with DataStore(f) as ds:
            return ds.sample_registry

    def _load_metadata(self):
        """Return metadata for first file.

        Assumes that metadata is identical for all feature files
        """
        metadata = dict()
        first_file = self.filenames[0]
        with DataStore(first_file) as ds:
            if ds._meta_group_ in ds.fh:
                for k in ds.fh[ds._meta_group_].keys():
                    metadata[k] = ds.get_meta(k)
        return metadata

    @property
    def index(self):
        """Return a dictionary describing all samples.

        The dictionary maps sample names to their file and HDF group. It is
        sorted by reference coordinate.
        """
        if self._index is None:
            self._index = self._get_sorted_index()
        return self._index

    def _get_sorted_index(self):
        """Get index of samples indexed by reference and ordered by start pos.

        :returns: {ref_name: [sample dicts sorted by start]}

        """
        ref_names = defaultdict(list)

        for sample_key, fn in self.samples:
            d = medaka.common.Sample.decode_sample_name(sample_key)
            if d is not None:
                d['sample_key'] = sample_key
                d['filename'] = fn
                ref_names[d['ref_name']].append(d)

        # sort dicts so that refs are in order and within a ref,
        # chunks are in order
        ref_names_ordered = OrderedDict()

        # sort by start and -end so that if we have two samples with the same
        # start but different end points, the longest sample comes first
        def get_major_minor(x):
            return tuple((int(i) for i in x.split('.')))

        def sorter(x):
            return get_major_minor(
                x['start']) + tuple((-i for i in get_major_minor(x['end'])))

        for ref_name in sorted(ref_names):
            ref_names[ref_name].sort(key=sorter)
            ref_names_ordered[ref_name] = ref_names[ref_name]

        return ref_names_ordered

    def yield_from_feature_files(self, regions=None, samples=None):
        """Yield `medaka.common.Sample` objects from one or more feature files.

        :regions: list of `medaka.common.Region` s for which to yield samples.
        :samples: iterable of sample names to yield (in order in which they
            are supplied).

        :yields: `medaka.common.Sample` objects.

        """
        if samples is not None:
            # yield samples in the order they are asked for
            for sample, fname in samples:
                yield DataStore(fname).load_sample(sample)
        else:
            all_samples = self.index
            if regions is None:
                regions = [
                    medaka.common.Region.from_string(x)
                    for x in sorted(all_samples)]
            for reg in regions:
                if reg.ref_name not in self.index:
                    continue
                for sample in self.index[reg.ref_name]:
                    # samples can have major.minor coords, round to end excl.
                    sam_reg = medaka.common.Region(
                        sample['ref_name'],
                        int(float(sample['start'])),
                        int(float(sample['end'])) + 1)
                    if sam_reg.overlaps(reg):
                        with DataStore(sample['filename']) as store:
                            yield store.load_sample(sample['sample_key'])
