"""Storing of training and inference data to file."""
from collections import defaultdict, OrderedDict
from concurrent.futures import \
    as_completed, ProcessPoolExecutor, ThreadPoolExecutor
import warnings

import dill
import numpy as np

import medaka.common

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class DataStore(object):
    """Read and write data to .hdf files."""

    # these metadata Datasets are common to data files storing
    # both sample and model information
    # they are placed in a self._meta_group_ Group
    _meta_group_ = 'meta'
    _metadata_datasets_ = (
        'feature_encoder',  # pickled FeatureEncoder object
        'model_function',   # pickled partial function; creates model object
        'label_scheme')     # pickled LabelScheme object

    _sample_path_ = 'samples/data'  # data group contains sample Datasets
    _sample_registry_path_ = 'samples/registry'  # set of sample keys

    def __init__(self, filename, mode='r', verify_on_close=True):
        """Initialize a datastore.

        :param filename: file to open.
        :param mode: file opening mode ('r', 'w', 'a').
        :param verify_on_close: on file close, check that all samples logged
            as being stored in file have a corresponding group within the
            `.hdf`."
        """
        self.filename = filename
        self.mode = mode
        self.verify_on_close = verify_on_close

        self.logger = medaka.common.get_named_logger('DataStore')

        self.write_executor = ThreadPoolExecutor(1)
        self.write_futures = []

    def __enter__(self):
        """Create filehandle."""
        self.fh = h5py.File(self.filename, self.mode)
        self._load_metadata()  # accessed via self.metadata
        self._load_sample_registry()  # accessed via self.sample_registry

        return self

    def __exit__(self, *args):
        """Verify file if requested."""
        if self.mode != 'r':
            if self.verify_on_close:
                self._verify()
            else:
                self.logger.debug("Skipping validation on close.")
            self._write_metadata()
            self._write_sample_registry()
            self.write_executor.shutdown(wait=True)
        self.fh.close()

    def _verify(self):
        """Remove samples from registry if nonexistent or missing fields."""
        self.logger.debug("Verifying data.")
        self.fh.flush()
        fh = h5py.File(self.filename, 'r')

        # ensure that sample registry only contains the keys of samples
        # that exist
        self._sync_sample_registry()

        # get union of all fields in all samples
        sample_fields = set()
        for key in self.sample_registry:
            sample_data_path = '/'.join((self._sample_path_, key))
            sample_fields.update(set(fh[sample_data_path]))

        # if a sample is missing fields, delete it from the registry
        for key in self.sample_registry:
            sample_data_path = '/'.join((self._sample_path_, key))
            missing_fields = sample_fields - set(fh[sample_data_path])
            if len(missing_fields):
                self.sample_registry.remove(key)
                self.logger.debug('Removing sample {} '.format(key) +
                                  'as {} not present.'.format(missing_fields))

    def write_sample(self, sample):
        """Write sample to hdf.

        Checks are performed to ensure a sample is not written twice and
        a count of unique training labels seen is maintained.

        :param sample: `medaka.common.Sample` object.
        """
        contains_numpy_array = any(isinstance(getattr(sample, field),
                                              np.ndarray)
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
            self.sample_registry.add(sample.name)
        else:
            self.logger.debug(
                'Not writing {} as present already'.format(
                    sample.name))

    def _write_dataset(self, location, data):
        """Write data, compressing numpy arrays."""
        if isinstance(data, np.ndarray):
            self.fh.create_dataset(location, data=data,
                                   compression='gzip',
                                   compression_opts=1)
        else:
            self.fh[location] = data

    def load_sample(self, key):
        """Load `medaka.common.Sample` object from file.

        :param key: str, sample name.
        :returns: `medaka.common.Sample` object.
        """
        s = dict()
        for field in medaka.common.Sample._fields:
            pth = '{}/{}/{}'.format(self._sample_path_, key, field)
            if pth in self.fh:
                s[field] = self.fh[pth][()]
                # handle loading of bytestrings
                if isinstance(s[field], np.ndarray) and \
                        isinstance(s[field][0], type(b'')):
                    s[field] = np.char.decode(s[field])
            else:
                s[field] = None
        return medaka.common.Sample(**s)

    def _load_pickled(self, path):
        """Load and return pickled object."""
        obj = dill.loads(self.fh[path][()])
        return obj

    def _write_pickled(self, obj, path):
        """Write a pickled object to file."""
        if path in self.fh:
            del self.fh[path]
        pickled_obj = np.string_(dill.dumps(obj))
        self.fh[path] = pickled_obj

    def _load_feature_encoder(self, path):
        """Load and return feature encoder."""
        obj = self._load_pickled(path)
        # set logger; pickle does not handle correctly
        obj.logger = medaka.common.get_named_logger('Feature')
        return obj

    @property
    def _metadata_loaders(self):
        """Return dict of metadata loaders."""
        loaders = {'feature_encoder': self._load_feature_encoder,
                   'model_function': self._load_pickled,
                   'label_scheme': self._load_pickled}
        return loaders

    @property
    def _metadata_writers(self):
        """Return dict of metadata writers."""
        writers = {'feature_encoder': self._write_pickled,
                   'model_function': self._write_pickled,
                   'label_scheme': self._write_pickled}
        return writers

    def _load_metadata_dataset(self, dataset):
        """Load dataset using appropriate loader."""
        try:
            loader = self._metadata_loaders[dataset]
        except KeyError:
            self.logger.debug('No metadata loader defined for {}'.format(
                dataset))
        else:
            path = '/'.join((self._meta_group_, dataset))
            dataset = loader(path)
            return dataset

    def _load_metadata(self, datasets=None):
        """Load meta data."""
        if datasets is None:
            datasets = self._metadata_datasets_
        metadata = dict()
        for d in datasets:
            try:
                metadata[d] = self._load_metadata_dataset(d)
            except Exception as e:
                self.logger.debug("Could not load {} from {}. {}.".format(
                    d, self.filename, e))
        self.metadata = metadata

    def _write_metadata(self):
        """Write meta data."""
        for dataset, data in self.metadata.items():
            self._write_metadata_dataset(dataset, data)

    def _write_metadata_dataset(self, dataset, data):
        """Write metadata."""
        self.logger.debug("Writing metadata.")
        try:
            writer = self._metadata_writers[dataset]
        except KeyError:
            self.logger.debug('No metadata writer defined for {}'.format(
                dataset))
        else:
            path = '/'.join((self._meta_group_, dataset))
            writer(data, path)

    def _load_sample_registry(self):
        """Load sample registry."""
        try:
            sample_keys = self._load_pickled(
                self._sample_registry_path_)
        except KeyError:
            sample_keys = set()
        finally:
            self.sample_registry = sample_keys

    def _write_sample_registry(self):
        """Write sample registry."""
        self.logger.debug("Writing sample registry.")
        if self._sample_registry_path_ in self.fh:
            del self.fh[self._sample_registry_path_]
        self._write_pickled(self.sample_registry,
                            self._sample_registry_path_)

    def _sync_sample_registry(self):
        """Find samples in file and update registry accordingly."""
        try:
            sample_keys = set(self.fh[self._sample_path_])
        except KeyError:
            self.logger.debug('No {} found in {}.'.format(
                self._sample_path_, self.filename))
            sample_keys = set()
        finally:
            self.sample_registry = sample_keys

    @property
    def n_samples(self):
        """Return the number of samples stored in file."""
        return len(self.sample_registry)


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
        first_file = self.filenames[0]
        with DataStore(first_file) as ds:
            return ds.metadata

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
        # start but differrent end points, the longest sample comes first
        def get_major_minor(x):
            return tuple((int(i) for i in x.split('.')))

        def sorter(x):
            return get_major_minor(
                x['start']) + tuple((-i for i in get_major_minor(x['end'])))

        for ref_name in sorted(ref_names):
            ref_names[ref_name].sort(key=sorter)
            ref_names_ordered[ref_name] = ref_names[ref_name]

        return ref_names_ordered

    def yield_from_feature_files(self, ref_names=None, samples=None):
        """Yield `medaka.common.Sample` objects from one or more feature files.

        :ref_names: iterable of str, only process these references.
        :samples: iterable of sample names to yield (in order in which they
            are supplied).

        :yields: `medaka.common.Sample` objects.

        """
        if samples is not None:
            # yield samples in the order they are asked for
            for sample, fname in samples:
                yield DataStore(fname).load_sample(sample)
        else:
            # yield samples sorted by ref_name and start
            if ref_names is None:
                ref_names = sorted(self.index)
            for ref_name in ref_names:
                for d in self.index[ref_name]:
                    with DataStore(d['filename']) as ds:
                        yield ds.load_sample(d['sample_key'])
