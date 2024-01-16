"""Storing of training and inference data to file."""
from abc import ABC, abstractmethod
from collections import defaultdict, OrderedDict
from concurrent.futures import \
    as_completed, ProcessPoolExecutor, ThreadPoolExecutor
import contextlib
import os
import pickle
import tarfile
import tempfile
import warnings

import numpy as np

import medaka.common


with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class BaseModelStore(ABC):
    """Base class for model store classes."""

    @abstractmethod
    def __init__(*args, **kwargs):
        """Initialize feature encoder."""
        raise NotImplementedError

    @abstractmethod
    def load_model(self, time_steps):
        """Load a model from hdf file/tensorflow directory."""
        raise NotImplementedError

    @abstractmethod
    def get_meta(self, key):
        """Retrieve a meta data item."""
        raise NotImplementedError

    @abstractmethod
    def copy_meta(self, key):
        """Copy meta data to hdf."""
        raise NotImplementedError


class ModelStore(BaseModelStore):
    """Read and write model and meta to a hdf file."""

    def __init__(self, filepath):
        """Initialize a Modelstore.

        :param filename: filepath to hdf file
        """
        self.filepath = filepath
        self.logger = medaka.common.get_named_logger('MdlStore')

    def __enter__(self):
        """Create context for handling a modelstore file."""
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """Exit context manager."""
        if exception_type is not None:
            self.logger.info('ModelStore exception {}'.format(exception_value))

    def load_model(self, time_steps=None):
        """Load a model from an .hdf file.

        :param time_steps: number of time points in RNN, `None` for dynamic.

        ..note:: keras' `load_model` cannot handle CuDNNGRU layers, hence this
            function builds the model then loads the weights.
        """
        with DataStore(self.filepath) as ds:
            self.logger.info('filepath {}'.format(self.filepath))
            model_partial_function = ds.get_meta('model_function')
            model = model_partial_function(time_steps=time_steps)
            model.load_weights(self.filepath)
        return model

    def get_meta(self, key):
        """Retrieve a meta data item.

        :param key: name of item to load.
        """
        with DataStore(self.filepath) as ds:
            return ds.get_meta(key)

    def copy_meta(self, other):
        """Copy meta data to hdf."""
        with DataStore(self.filepath) as ds:
            return ds.copy_meta(other)


class ModelStoreTF(BaseModelStore):
    """Read and write model to tensorflow storage directory."""

    top_level_dir = 'model'

    def __init__(self, filepath):
        """Initialize a Modelstore.

        :param filename: filepath to saved_model directory
        """
        self.logger = medaka.common.get_named_logger('MdlStrTF')
        self.filepath = filepath
        self.meta = None
        self.tmpdir = None

    def unpack(self):
        """Unpack model files from archive."""
        if self.tmpdir is None:
            # tmpdir is removed by .cleanup()
            self.tmpdir = tempfile.TemporaryDirectory()
            self._exitstack = contextlib.ExitStack()
            self._exitstack.enter_context(self.tmpdir)
            with tarfile.open(self.filepath) as tar:
                # CVE-2007-4559. We only really extract our own files, but
                # someone could unwittingly run a model from the interwebs
                # and blame us when things go wrong.

                def _is_within_directory(directory, target):
                    abs_directory = os.path.abspath(directory)
                    abs_target = os.path.abspath(target)
                    prefix = os.path.commonprefix([abs_directory, abs_target])
                    return prefix == abs_directory

                def _safe_extract(
                        tar, path=".", members=None, *, numeric_owner=False):
                    for member in tar.getmembers():
                        member_path = os.path.join(path, member.name)
                        if not _is_within_directory(path, member_path):
                            raise Exception(
                                "Attempted Path Traversal in Tar File")
                    tar.extractall(path, members, numeric_owner=numeric_owner)

                # Do the extraction
                _safe_extract(tar, path=self.tmpdir.name)

            meta_file = os.path.join(
                self.tmpdir.name, self.top_level_dir, 'meta.pkl')
            with open(meta_file, 'rb') as fh:
                self.meta = pickle.load(fh)
        return self

    def cleanup(self):
        """Clean up temporary files."""
        if self.tmpdir:
            try:
                self._exitstack.close()
            except OSError:
                self.logger.warning(
                    "Failed to correctly clean up model temporary files. "
                    f"Some files might be left in {self.tmpdir.name}.")
            else:
                self.logger.info(
                    "Successfully removed temporary files "
                    f"from {self.tmpdir.name}.")
            self.meta = None
            self.tmpdir = None
            del self._exitstack

    def __enter__(self):
        """Context manager."""
        self.unpack()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """Remove temporary unpack_filepath."""
        self.cleanup()
        if exception_type is not None:
            self.logger.info('ModelStoreTF exception {}'.format(
                exception_type))

    def __del__(self):
        """Run cleanup on destroy."""
        self.cleanup()

    def load_model(self, time_steps=None):
        """Load a model from a tf saved_model file.

        :param time_steps: number of time points in RNN, `None` for dynamic.

        ..note:: this function builds the model then loads the weights.
        """
        self.unpack()
        model_partial_function = self.get_meta('model_function')
        self.model = model_partial_function(time_steps=time_steps)
        self.logger.info("Model {}".format(self.model))
        weights = os.path.join(
            self.tmpdir.name, self.top_level_dir, 'variables', 'variables')
        self.logger.info(
            "loading weights from {} (using expect partial)".format(weights))
        # expect partial ignores errors about the optimizer state not being
        # present saving the optimizer state would make the models bigger.
        # would be nice to figure out how to delete the references to the
        # optimizer from the models
        self.model.load_weights(weights).expect_partial()
        return self.model

    def get_meta(self, key):
        """Load (deserialise) a meta data item.

        :param key: name of item to load.
        """
        self.unpack()
        return self.meta[key]

    def copy_meta(self, hdf):
        """Copy metadata to hdf file.

        :param hdf: filename of hdf file.
        """
        self.unpack()
        with DataStore(hdf, 'a') as ds:
            for k, v in self.meta.items():
                ds.set_meta(v, k)


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

        self.logger = medaka.common.get_named_logger('DataStre')

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
                            isinstance(data[0], np.compat.unicode):
                        data = np.char.encode(data)
                    location = '{}/{}/{}'.format(
                        self._sample_path_, sample.name, field)
                    self.write_futures.append(
                        self.write_executor.submit(
                            self._write_dataset, location, data))
            self.logger.debug(
                'Adding {} to sample registry'.format(
                    sample.name))
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
        s = {x: None for x in medaka.common.Sample._fields}
        group = self.fh['{}/{}'.format(self._sample_path_, key)]
        for field in medaka.common.Sample._fields:
            try:
                s[field] = group[field][()]
            except KeyError:
                pass
            else:
                # handle loading of bytestrings
                if isinstance(s[field], np.ndarray) and \
                        isinstance(s[field][0], type(b'')):
                    s[field] = np.char.decode(s[field])
                if isinstance(s[field], bytes):
                    s[field] = s[field].decode()
        return medaka.common.Sample(**s)

    def _write_dataset(self, location, data):
        """Write data, compressing numpy arrays."""
        if isinstance(data, np.ndarray):
            self.fh.create_dataset(
                location, data=data, compression='gzip', compression_opts=1)
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
                self.logger.debug("Loaded sample register.")
            except KeyError:
                if self._sample_path_ in self.fh:
                    sample_registry = set(self.fh[self._sample_path_].keys())
                    if sample_registry:
                        self._sample_registry = sample_registry
                        self.logger.debug("Created missing sample register.")
                else:
                    self._sample_registry = set()

    def _write_sample_registry(self):
        """Write sample registry."""
        self.logger.debug("Writing sample registry.")
        if self._sample_registry_path_ in self.fh:
            del self.fh[self._sample_registry_path_]
        self._write_pickled(
            self.sample_registry, self._sample_registry_path_)


class DataIndex(object):
    """Index and serve samples from multiple `DataStore` compatible files."""

    def __init__(self, filenames, threads=4):
        """Intialize an index across a set of files.

        :param filenames: list of files to index.
        :param threads: number of threads to use for indexing.
        """
        self.logger = medaka.common.get_named_logger('DataIndx')
        if isinstance(filenames, str):
            filenames = [filenames]
        self.filenames = filenames
        self.threads = threads
        self.n_files = len(self.filenames)
        self._index = None
        self._extract_sample_registries()
        self.metadata = self._load_metadata()
        # in cases where where we have a single file containing samples over
        # many contigs we might as well open the file once rather than once
        # per contig.
        self._ds = DataStore(filenames[0])

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

    @property
    def regions(self):
        """Return a Region object for each ref_name in the index.

        Each region will have start and end set to None.
        """
        return [medaka.common.Region(r, None, None) for r in self.index]

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
        if samples is None:
            all_samples = self.index
            if regions is None:
                regions = [
                    medaka.common.Region.from_string(x)
                    for x in sorted(all_samples)]
            samples = list()
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
                        samples.append(
                            (sample['sample_key'], sample['filename']))
        # yield samples reusing filehandle where possible
        for key, fname in samples:
            if fname != self._ds.filename:
                self._ds = DataStore(fname)
            yield self._ds.load_sample(key)
