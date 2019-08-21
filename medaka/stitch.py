"""Creation of contiguous consensus sequences from chunked network outputs."""

import itertools

import medaka.common
import medaka.datastore
import medaka.labels


def write_fasta(filename, contigs):
    """Write a fasta file from tuples of (name, sequence).

    :param filename: output filename.
    :param contigs: tuples of the form (sequence name, base sequence).

    """
    with open(filename, 'w') as fasta:
        for name, seq in contigs:
            fasta.write('>{}\n{}\n'.format(name, seq))


def stitch_from_probs(h5_fp, regions=None):
    """Join overlapping label probabilities from HDF5 files.

     Network outputs from multiple samples stored within a file are spliced
     together into a logically contiguous array and decoded to generate
     contiguous sequence(s).

    :param h5_fp: iterable of HDF5 filepaths
    :param regions: iterable of regions to process

    :returns: list of (region string, sequence)
    """
    logger = medaka.common.get_named_logger('Stitch')

    index = medaka.datastore.DataIndex(h5_fp)

    # lookup LabelScheme name stored in HDF5 as piece of metadata
    try:
        label_class = medaka.labels.label_schemes[
            index.meta['label_scheme_class']]
    except KeyError:
        logger.debug(
            "Could not find `label_scheme_class` definition in input file, "
            "assuming HaploidLabelScheme.")
        label_class = medaka.labels.HaploidLabelScheme
    finally:
        label_scheme = label_class()

    logger.debug("Label decoding is:\n{}".format(
        '\n'.join('{}: {}'.format(k, v)
                  for k, v in label_scheme._decoding.items())))

    if regions is None:
        ref_names = index.index.keys()
    else:
        # TODO: respect entire region specification
        ref_names = list()
        for r in regions:
            region = medaka.common.Region.from_string(r)
            if region.start is None or region.end is None:
                logger.warning("Ignoring start:end for '{}'.".format(region))
            ref_names.append(region.ref_name)

    def get_pos(sample, i):
        return '{}.{}'.format(
            sample.positions[i]['major'] + 1, sample.positions[i]['minor'])

    ref_assemblies = []
    for ref_name in ref_names:
        logger.info("Processing {}.".format(ref_name))
        data_gen = index.yield_from_feature_files(ref_names=(ref_name,))
        seq_parts = list()
        # first sample
        s1 = next(data_gen)
        start = get_pos(s1, 0)
        start_1 = None
        start_2 = None

        for s2 in itertools.chain(data_gen, (None,)):
            s1_name = 'Unknown' if s1 is None else s1.name
            s2_name = 'Unknown' if s2 is None else s2.name

            # s1 is last chunk
            if s2 is None:
                end_1 = None
            else:
                # s2 ends before s1
                if s2.last_pos <= s1.last_pos:
                    logger.info('{} ends before {}, skipping.'.format(
                        s2_name, s1_name
                    ))
                    continue
                # s1 and s2 overlap by only one position
                # or there is no overlap between s1 and s2
                elif s2.first_pos >= s1.last_pos:
                    # trigger a break
                    end_1, start_2 = None, None
                else:
                    try:
                        end_1, start_2 = medaka.common.Sample.overlap_indices(
                            s1, s2)
                    except medaka.common.OverlapException as e:
                        logger.info(
                            "Unhandled overlap type whilst stitching chunks.")
                        raise(e)

            new_seq = label_scheme.decode_consensus(
                s1.slice(slice(start_1, end_1)))

            seq_parts.append(new_seq)

            if end_1 is None:
                region_string = '{}:{}-{}'.format(s1.ref_name, start,
                                                  get_pos(s1, -1))
                ref_assemblies.append((region_string, ''.join(seq_parts)))
                seq_parts = list()

                if s2 is not None and start_2 is None:
                    msg = 'There is no overlap betwen {} and {}'
                    logger.info(msg.format(s1_name, s2_name))
                    start = get_pos(s2, 0)

            s1 = s2
            start_1 = start_2
    return ref_assemblies


def stitch(args):
    """Entry point for stitching program."""
    joined = stitch_from_probs(args.inputs, regions=args.regions)
    write_fasta(args.output, joined)
