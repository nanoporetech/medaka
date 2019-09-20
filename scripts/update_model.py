#!/usr/bin/env python3
import argparse
import functools
import os
import pickle
import shutil
import sys
import yaml

import h5py
import numpy as np

import medaka.features
import medaka.labels
import medaka.models


def main():

    parser = argparse.ArgumentParser(
        description='Update a pre-version-0.9.0 model.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('model', help='input model filepath')
    parser.add_argument('output', help='output model filepath')
    args = parser.parse_args()

    if os.path.exists(args.output):
        sys.stderr.write('{} exists already!\n'.format(args.output))
        sys.exit(1)

    shutil.copy(args.model, args.output)

    with h5py.File(args.output) as h:

        # check that model can be built
        model_name = yaml.unsafe_load(h['medaka_model_name'][()])
        try:
            build_model = medaka.models.model_builders[model_name]
        except KeyError('Can not convert; requires deprecated ' + \
                        '{} function.'.format(model_name)):
            sys.exit(1)

        features = yaml.unsafe_load(h['medaka_feature_decoding'][()])
        feat_len = len(features)
        classes = yaml.unsafe_load(h['medaka_label_decoding'][()])
        num_classes = len(classes)

        gru_size = 128
        classify_activation = 'softmax'
        # load specified model kwargs if they exist
        model_kwargs = yaml.unsafe_load(h['medaka_model_kwargs'][()])
        if 'gru_size' in model_kwargs:
            gru_size = model_kwargs['gru_size']
        if 'classify_activation' in model_kwargs:
            activation = model_kwargs['classify_activation']

        normalise = 'total'
        medaka_features_kwargs = yaml.unsafe_load(
            h['medaka_features_kwargs'][()])
        if 'normalise' in medaka_features_kwargs:
            normalise = medaka_features_kwargs['normalise']

        # delete existing metadata
        for i in ['medaka_feature_decoding',
                  'medaka_features_kwargs',
                  'medaka_label_counts',
                  'medaka_label_decoding',
                  'medaka_model_kwargs',
                  'medaka_model_name']:
            if h.get(i):
                del h[i]

    # write new-style metadata
    with medaka.datastore.DataStore(args.output, mode='a') as ds:

        ds.set_meta(medaka.labels.HaploidLabelScheme(), 'label_scheme')
        ds.set_meta(
            medaka.features.CountsFeatureEncoder(normalise=normalise),
            'feature_encoder')
        ds.set_meta(
            functools.partial(
                build_model, feat_len, num_classes,
                gru_size=gru_size, classify_activation=classify_activation),
            'model_function')

if __name__ == '__main__':
    main()

