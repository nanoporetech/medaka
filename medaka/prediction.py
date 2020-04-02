"""Inference program and ancilliary functions."""
from concurrent.futures import as_completed, ThreadPoolExecutor
import logging
import os
import queue
import threading
from timeit import default_timer as now

import numpy as np

import medaka.common
import medaka.datastore
import medaka.features
import medaka.models


def run_prediction(
        output, bam, regions, model, feature_encoder,
        chunk_len, chunk_ovlp, batch_size=200,
        save_features=False, enable_chunking=True):
    """Inference worker."""
    logger = medaka.common.get_named_logger('PWorker')

    remainder_regions = list()
    loader = DataLoader(
        4 * batch_size, bam, regions, feature_encoder,
        chunk_len=chunk_len, chunk_overlap=chunk_ovlp,
        enable_chunking=enable_chunking)
    batches = medaka.common.grouper(loader, batch_size)

    total_region_mbases = sum(r.size for r in regions) / 1e6
    logger.info(
        "Running inference for {:.1f}M draft bases.".format(
            total_region_mbases))

    with medaka.datastore.DataStore(output, 'a') as ds:
        mbases_done = 0
        cache_size_log_interval = 5

        t0 = now()
        tlast = t0
        tcache = t0
        for data in batches:
            if now() - tcache > cache_size_log_interval:
                logger.info("Samples in cache: {}.".format(
                    loader.results.qsize()))
                tcache = now()
            x_data = np.stack([x.features for x in data])
            class_probs = model.predict_on_batch(x_data)
            # calculate bases done taking into account overlap
            new_bases = 0
            for x in data:
                if chunk_ovlp < x.size:
                    new_bases += x.last_pos[0] - x._get_pos(chunk_ovlp)[0]
                else:
                    new_bases += x.span
            mbases_done += new_bases / 1e6
            # just to avoid funny log msg...
            mbases_done = min(mbases_done, total_region_mbases)
            t1 = now()
            if t1 - tlast > 10:
                tlast = t1
                msg = '{:.1%} Done ({:.1f}/{:.1f} Mbases) in {:.1f}s'
                logger.info(msg.format(
                    mbases_done / total_region_mbases, mbases_done,
                    total_region_mbases, t1 - t0))

            for sample, prob, feat in zip(data, class_probs, x_data):
                # write out positions and predictions for later analysis
                features = feat if save_features else None
                new_sample = sample.amend(
                    label_probs=prob, features=features)
                ds.write_sample(new_sample)

    remainder_regions = loader.remainders
    logger.info("All done, {} remainder regions.".format(
        len(remainder_regions)))
    return remainder_regions


def predict(args):
    """Inference program."""
    logger_level = logging.getLogger(__package__).level
    if logger_level > logging.DEBUG:
        os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

    import tensorflow as tf
    from tensorflow.keras import backend as K

    args.regions = medaka.common.get_regions(
        args.bam, region_strs=args.regions)
    logger = medaka.common.get_named_logger('Predict')
    logger.info('Processing region(s): {}'.format(
        ' '.join(str(r) for r in args.regions)))

    # create output and copy meta
    with medaka.datastore.DataStore(args.model) as ds:
        ds.copy_meta(args.output)
        feature_encoder = ds.get_meta('feature_encoder')

    feature_encoder.tag_name = args.tag_name
    feature_encoder.tag_value = args.tag_value
    feature_encoder.tag_keep_missing = args.tag_keep_missing
    feature_encoder.read_group = args.RG

    logger.info("Setting tensorflow threads to {}.".format(args.threads))
    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
    K.set_session(tf.Session(
        config=tf.ConfigProto(
            intra_op_parallelism_threads=args.threads,
            inter_op_parallelism_threads=args.threads)
    ))
    if tf.test.is_gpu_available(cuda_only=True):
        logger.info("Found a GPU.")
        logger.info(
            "If cuDNN errors are observed, try setting the environment "
            "variable `TF_FORCE_GPU_ALLOW_GROWTH=true`. To explicitely "
            "disable use of cuDNN use the commandline option "
            "`--disable_cudnn. If OOM (out of memory) errors are found "
            "please reduce batch size.")

    # Split overly long regions to maximum size so as to not create
    #   massive feature matrices
    MAX_REGION_SIZE = int(1e6)  # 1Mb
    regions = []
    for region in args.regions:
        if region.size > MAX_REGION_SIZE:
            # chunk_ovlp is mostly used in overlapping pileups (which generally
            # end up being expanded compared to the draft coordinate system)
            regs = region.split(
                MAX_REGION_SIZE, overlap=args.chunk_ovlp, fixed_size=False)
        else:
            regs = [region]
        regions.extend(regs)

    logger.info("Processing {} long region(s) with batching.".format(
        len(regions)))

    logger.info("Using model: {}.".format(args.model))

    model = medaka.models.load_model(args.model, time_steps=args.chunk_len,
                                     allow_cudnn=args.allow_cudnn)

    # the returned regions are those where the pileup width is smaller than
    # chunk_len
    remainder_regions = run_prediction(
        args.output, args.bam, regions, model, feature_encoder,
        args.chunk_len, args.chunk_ovlp,
        batch_size=args.batch_size, save_features=args.save_features
    )

    # short/remainder regions: just do things without chunking. We can do this
    # here because we now have the size of all pileups (and know they are
    # small).
    # TODO: can we avoid calculating pileups twice whilst controlling memory?
    if len(remainder_regions) > 0:
        logger.info("Processing {} short region(s).".format(
            len(remainder_regions)))
        model = medaka.models.load_model(args.model, time_steps=None,
                                         allow_cudnn=args.allow_cudnn)
        for region in remainder_regions:
            new_remainders = run_prediction(
                args.output, args.bam, [region[0]], model, feature_encoder,
                args.chunk_len, args.chunk_ovlp,  # these won't be used
                batch_size=args.batch_size, save_features=args.save_features,
                enable_chunking=False
            )
            if len(new_remainders) > 0:
                # shouldn't get here
                ignored = [x[0] for x in new_remainders]
                n_ignored = len(ignored)
                logger.warning("{} regions were not processed: {}.".format(
                    n_ignored, ignored))

    logger.info("Finished processing all regions.")

    if args.check_output:
        logger.info("Validating and finalising output data.")
        with medaka.datastore.DataStore(args.output, 'a') as ds:
            pass


class DataLoader(object):
    """Loading of data for inference."""

    def __init__(self, sample_cache_size, bam, regions, *args, **kwargs):
        """Initialise data loading.

        :param sample_cache_size: maximum number of network inputs to
            pre-load.
        :param bam: input `.bam` file.
        :param regions: regions to process.
        :param args: position arguments to use when creating
            `features.SampleGenerator` instances.
        :param kwargs: keyword arguments to use when creating
            `features.SampleGenerator` instances.
        """
        self.logger = medaka.common.get_named_logger('DLoader')
        self.sample_cache_size = sample_cache_size
        self.region_cache_size = 4
        self.bam = bam
        self.regions = regions
        self.args = args
        self.kwargs = kwargs
        self.results = queue.Queue(maxsize=sample_cache_size)
        self.have_data = threading.Event()
        self.remainders = []

        self.logger.info('Initializing data loader')
        self.have_data.set()
        self.thread = threading.Thread(target=self._fill_parallel)
        self.thread.daemon = True
        self.thread.start()

    def _fill_serial(self):
        # process one region at a time
        for region in self.regions:
            samples, remain = self._run_region(
                self.bam, region, *self.args, **self.kwargs)
            for sample in samples:
                self.results.put(sample)
            self.remainders.extend(remain)
        self.have_data.clear()

    def _fill_parallel(self):
        # process multiple regions at a time, up to a maximum to limit memory
        # use. Note that the number of workers also serves as a memory
        # limit when data is being consumed as fast as it is produced.
        regions = iter(self.regions)
        futures = dict()
        submitted = True
        t0 = now()
        cache_check_interval = 3
        min_region_cache = 4
        with ThreadPoolExecutor(max_workers=4) as executor:
            while True:
                if submitted:
                    try:
                        submit_reg = next(regions)
                    except StopIteration:
                        break
                # try to submit
                if len(futures) < self.region_cache_size:
                    self.logger.debug("Submitting {}.".format(submit_reg))
                    futures[str(submit_reg)] = executor.submit(
                        self._run_region, self.bam, submit_reg,
                        *self.args, **self.kwargs)
                    submitted = True
                else:
                    submitted = False
                # try to fetch
                done = []
                for kreg, fut in futures.items():
                    if fut.done():
                        samples, remain = fut.result()
                        self.remainders.extend(remain)
                        for sample in samples:
                            self.results.put(sample)
                        done.append(kreg)
                for kreg in done:
                    del futures[kreg]
                # keep things flowing
                if now() - t0 > cache_check_interval:
                    t0 = now()
                    if self.results.qsize() < 0.5 * self.sample_cache_size:
                        self.logger.debug(
                            "Expanding region cache from {},".format(
                                self.region_cache_size))
                        self.region_cache_size += 1
                    elif self.results.qsize() > 0.9 * self.sample_cache_size:
                        self.logger.debug(
                            "Reducing region cache from {},".format(
                                self.region_cache_size))
                        self.region_cache_size = max(
                            min_region_cache, self.region_cache_size - 1)
                    else:
                        self.logger.debug("Region cache is good size.")
            # collect remaining futures
            for fut in as_completed(futures.values()):
                samples, remain = fut.result()
                self.remainders.extend(remain)
                for sample in samples:
                    self.results.put(sample)
        # signal everything has been processed
        self.have_data.clear()

    @staticmethod
    def _run_region(bam, region, *args, **kwargs):
        data_gen = medaka.features.SampleGenerator(
            bam, region, *args, **kwargs)
        return data_gen.samples, data_gen._quarantined

    def __iter__(self):
        """Simply return self."""
        return self

    def __next__(self):
        """Iterate over items in internal network input cache."""
        while self.have_data.is_set():
            try:
                res = self.results.get(timeout=0.1)
            except queue.Empty:
                if not self.have_data.is_set():
                    break
            else:
                return res

        if not self.results.empty():
            return self.results.get()
        raise StopIteration
