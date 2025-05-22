"""Inference program and ancilliary functions."""
import queue
import threading
from timeit import default_timer as now

from medaka.architectures.base_classes import ReadLevelFeaturesModel
import medaka.common
import medaka.datastore
import medaka.features
import medaka.models
import medaka.torch_ext


def run_prediction(
        output, bam, regions, model, feature_encoder,
        chunk_len, chunk_ovlp, batch_size=200,
        save_features=False, enable_chunking=True,
        bam_workers=2):
    """Inference worker."""
    logger = medaka.common.get_named_logger('PWorker')

    # start loading data
    loader = DataLoader(
        bam, regions, batch_size,
        batch_cache_size=8, bam_workers=bam_workers,
        feature_encoder=feature_encoder,
        chunk_len=chunk_len, chunk_overlap=chunk_ovlp,
        enable_chunking=enable_chunking)

    total_region_mbases = sum(r.size for r in regions) / 1e6
    logger.info(
        "Running inference for {:.1f}M draft bases.".format(
            total_region_mbases))

    remainder_regions = list()
    n_batches = 0
    with medaka.datastore.DataStore(output, 'a') as ds:
        mbases_done = 0
        cache_size_log_interval = 5

        t0 = now()
        tlast = t0
        tcache = t0
        for data, batch in loader:
            n_batches += 1
            class_probs = model.predict_on_batch(batch)
            for sample, prob, feat in zip(data, class_probs, batch.features):
                # write out positions and predictions for later analysis
                features = feat if save_features else None
                new_sample = sample.amend(
                    label_probs=prob, features=features)
                ds.write_sample(new_sample)

            # log loading of batches
            if now() - tcache > cache_size_log_interval:
                logger.debug("Batches/Samples in cache: {}/{}".format(
                    loader._batches.qsize(), loader._samples.qsize()))
                tcache = now()
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
                logger.debug(msg.format(
                    mbases_done / total_region_mbases, mbases_done,
                    total_region_mbases, t1 - t0))

    remainder_regions = loader.remainders
    logger.info("Processed {} batches".format(n_batches))
    logger.info("All done, {} remainder regions.".format(
        len(remainder_regions)))
    return remainder_regions


def predict(args):
    """Inference program."""
    logger = medaka.common.get_named_logger('Predict')

    bam_regions = medaka.common.get_bam_regions(
        args.bam, regions=args.regions)
    logger.info('Processing region(s): {}'.format(
        ' '.join(str(r) for r in bam_regions)))

    # split out regions which are smaller than the chunk
    # size for processing later.
    regions, remainder_regions = [], []
    for region in bam_regions:
        if region.size < args.chunk_len:
            remainder_regions.append(region)
        else:
            # Split overly long regions to maximum size so as to not create
            #   massive feature matrices
            if region.size > args.bam_chunk:
                # chunk_ovlp is mostly used in overlapping pileups (which
                # generally end up being expanded compared to the draft
                # coordinate system)
                regions.extend(region.split(
                    args.bam_chunk, overlap=args.chunk_ovlp,
                    fixed_size=False))
            else:
                regions.append(region)

    logger.info("Using model: {}.".format(args.model))
    with medaka.models.open_model(args.model) as model_store:
        feature_encoder = model_store.get_meta('feature_encoder')
        model_store.copy_meta(args.output)

        feature_encoder.tag_name = args.tag_name
        feature_encoder.tag_value = args.tag_value
        feature_encoder.tag_keep_missing = args.tag_keep_missing
        feature_encoder.read_group = args.RG
        # If min_mapq was not set in the feature encoder, set it to the old
        # default
        if not hasattr(feature_encoder, "min_mapq"):
            feature_encoder.min_mapq = 1
        # Override feature_encoder min_mapq if it has been explicitly set
        if args.min_mapq is not None:
            feature_encoder.min_mapq = args.min_mapq
        msg = "Using minimum mapQ threshold of {} for read filtering."
        logger.info(msg.format(feature_encoder.min_mapq))

        # To support models legacy
        if not hasattr(feature_encoder, "sym_indels"):
            feature_encoder.sym_indels = False

        import torch
        if torch.cuda.device_count() > 0 and not args.cpu:
            logger.info("Found a GPU.")
            device = torch.device("cuda")
            # logger.info(
            #     "If cuDNN errors are observed, try setting the environment "
            #     "variable `TF_FORCE_GPU_ALLOW_GROWTH=true`. To explicitely "
            #     "disable use of cuDNN use the commandline option "
            #     "`--disable_cudnn. If OOM (out of memory) errors are found "
            #     "please reduce batch size.")
        else:
            device = torch.device("cpu")
            torch.set_num_threads(args.threads)
            args.full_precision = True

        model = model_store.load_model(device=device)

        # fails if invalid feature encoder for model
        model.check_feature_encoder_compatibility(feature_encoder)

        if isinstance(model, ReadLevelFeaturesModel) and device.type == "cpu":
            msg = (
                "WARNING: Using a read level features model on a CPU. This is "
                "expected to be slow. To improve performance, consider "
                "running on a GPU or using a counts matrix model, by "
                "passing the '--fast' to medaka_consensus.")
            logger.warning(msg)

    logger.info("Model device: {}".format(model.device()))
    if args.full_precision:
        logger.info("Running prediction at full precision")
    else:
        logger.info("Running prediction at half precision")
        model.half()

    # # TODO: see whether compilation improves performance
    # compiled_model = torch.compile(model, backend="inductor")

    bam_pool = medaka.features.BAMHandler(args.bam)

    if len(regions) > 0:

        logger.info("Processing {} long region(s) with batching.".format(
            len(regions)))

        # the returned regions are those where the pileup width is smaller
        # than chunk_len (which could happen due to gaps in coverage)
        remainder_regions_depth = run_prediction(
            args.output, bam_pool, regions, model, feature_encoder,
            args.chunk_len, args.chunk_ovlp,
            batch_size=args.batch_size, save_features=args.save_features,
            bam_workers=args.bam_workers)

        # run_prediction returns [(region, pileup width)]
        remainder_regions.extend([r[0] for r in remainder_regions_depth])

    # short/remainder regions: just do things without chunking. We can do
    # this here because we now have the size of all pileups (and know they
    # are small).
    # TODO: can we avoid calculating pileups twice whilst controlling
    # memory?
    if len(remainder_regions) > 0:
        logger.info("Processing {} short region(s).".format(
            len(remainder_regions)))

        # TODO: see whether this is still relevant for pytorch models
        # if "keras" in model.__module__:
        #     model.run_eagerly = True

        new_remainders = run_prediction(
            args.output, bam_pool, remainder_regions, model,
            feature_encoder,
            args.chunk_len, args.chunk_ovlp,  # these won't be used
            batch_size=1,  # everything is a different size, cant batch
            save_features=args.save_features, enable_chunking=False)
        if len(new_remainders) > 0:
            # shouldn't get here
            ignored = [x[0] for x in new_remainders]
            n_ignored = len(ignored)
            logger.warning("{} regions were not processed: {}.".format(
                n_ignored, ignored))

    logger.info("Finished processing all regions.")

    if args.check_output:
        logger.info("Validating and finalising output data.")
        with medaka.datastore.DataStore(args.output, 'a') as _:
            pass


class DataLoader(object):
    """Loading of data for inference."""

    def __init__(
            self, bam, regions, batch_size, batch_cache_size=8,
            bam_workers=2, **kwargs):
        """Initialise data loading.

        Once constructed, iterating over this object will yield
        tuples containing a list of `Samples`, and a corresponding
        medaka.torch_ext.Batch

        :param bam: input `.bam` file.
        :param regions: regions to process.
        :param batch_size: number samples in an inference batch.
        :param batch_cache_size: maximum number of inference batches
            to preload.
        :param bam_workers: number of worker threads for `.bam` reading.
        :param kwargs: keyword arguments to use when creating
            `features.SampleGenerator` instances.
        """
        self.logger = medaka.common.get_named_logger('DLoader')
        self.bam = bam
        self.batch_size = batch_size
        self.bam_workers = bam_workers
        self.kwargs = kwargs

        # Memory use is dominated by number of workers as each will
        # produce many batches in one shot. We make workers
        # block when many samples (or batches) have been calculated.
        # This setup will lead to at most 2 * batch_cache_size being
        # held in the system, plus anyhing computed but not transferred
        # into the self._samples
        sample_cache_size = batch_cache_size * batch_size
        self._samples = queue.Queue(maxsize=sample_cache_size)
        self._batches = queue.Queue(maxsize=batch_cache_size)
        # fill input queue before workers start
        self._regions = queue.Queue()
        for r in regions:
            self._regions.put(r)
        self.remainders = list()

        # loading of feature data into samples
        self._sample_workers = list()
        for i in range(bam_workers):
            t = threading.Thread(
                target=self._region_worker, name="Sample-{}".format(i))
            t.daemon = True
            t.start()
            self._sample_workers.append(t)
        # loading of samples into batches
        self._bthread = threading.Thread(
            target=self._batch_worker, name="Batcher")
        self._bthread.daemon = True
        self._bthread.start()

        self.logger.debug(
            "Initialised. Batch size:{}. Workers:{}. Cache:{}".format(
                batch_size, bam_workers, sample_cache_size))

    def __iter__(self):
        """Simply return self."""
        return self

    def __next__(self):
        """Generate batches of data for inference.

        :returns: (list of `Samples`, ndarray of stacks features).
        """
        while True:
            try:
                batch = self._batches.get_nowait()
                self.logger.debug(
                    "Batches ready: {}. Samples ready: {} ({} batches)".format(
                        self._batches.qsize(), self._samples.qsize(),
                        self._samples.qsize()/self.batch_size))
            except queue.Empty:
                pass
            else:
                self._batches.task_done()
                if batch is StopIteration:
                    # we have a single batch worker, safe to stop after
                    # seeing a single StopIteration
                    break
                return batch
        raise StopIteration

    @staticmethod
    def _run_region(bam, region, *args, **kwargs):
        data_gen = medaka.features.SampleGenerator(
            bam, region, *args, **kwargs)
        return data_gen.samples, data_gen._quarantined

    def _region_worker(self):
        t = threading.current_thread()
        self.logger.debug("{}: started".format(t.name))
        while True:
            try:
                region = self._regions.get_nowait()
            except queue.Empty:
                # when queue is empty there will never be more items
                self._samples.put(StopIteration)
                break
            else:
                samples, remain = self._run_region(
                    self.bam, region, **self.kwargs)
                for sample in samples:
                    self._samples.put(sample)
                self.remainders.extend(remain)
                self._regions.task_done()
        self.logger.debug("{}: finished".format(t.name))

    def _get_samples(self):
        """Iterate over items in internal network input cache."""
        nstops = 0  # careful to track workers have all signalled stop
        while True:
            try:
                res = self._samples.get_nowait()
            except queue.Empty:
                # more items may come
                pass
            else:
                if res is StopIteration:
                    nstops += 1
                    if nstops == self.bam_workers:
                        self._samples.task_done()
                        break
                else:
                    yield res
                    self._samples.task_done()

    def _batch_worker(self):
        """Create batch from list of samples to be used for inference.

        :returns: (list of `Samples`, medaka.torch_ext Batch).

        The first element of output used to keep track of which sample comes
        from where in the assembly, the second element is the actual data
        passed to the model's predict_on_batch function.
        """
        for data in medaka.common.grouper(
                self._get_samples(), self.batch_size):

            batch = medaka.torch_ext.Batch.collate(data)
            self._batches.put((data, batch))
        self._batches.put(StopIteration)
