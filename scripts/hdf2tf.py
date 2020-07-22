import argparse
import os
import pickle
import sys

from tensorflow.keras.models import save_model

import medaka.models


def main():

    parser = argparse.ArgumentParser(
        description='Convert hdf5 model file to tensorflow model',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_model', help='Input model filepath')
    parser.add_argument('output_model', 
                        help='Output model filepath without .tar.gz extension') 
    parser.add_argument('--force', action='store_true', 
                        help='Allow existing model to be overwritten')
    args = parser.parse_args()

    if os.path.exists(args.output_model + '.tar.gz') and  not args.force:
        sys.stderr.write('{} exists already. Use --force to overwrite.\n'
                    .format(args.output_model))
        sys.exit(1)

    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

    with medaka.datastore.ModelStore(args.input_model) as ms:
        partial_model_function = ms.get_meta('model_function')
        feature_encoder = ms.get_meta('feature_encoder')
        label_scheme = ms.get_meta("label_scheme")
        model = ms.load_model(time_steps=None)

    medaka_meta = {'model_function': partial_model_function,
                  'label_scheme': label_scheme,
                  'feature_encoder': feature_encoder}

    save_model(model, args.output_model)
    meta_filepath = os.path.join(args.output_model, 'meta.pkl')
    with open(meta_filepath, 'wb') as handle:
        pickle.dump(medaka_meta, handle)

    medaka.datastore.tar_dir(args.output_model, args.output_model + '.tar.gz')
    sys.stdout.write("Model written to {}\n".format(args.output_model + '.tar.gz'))
    medaka.datastore.del_dir(args.output_model)


if __name__ == '__main__':
    main()

