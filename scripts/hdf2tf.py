import argparse
import os
import pickle
import shutil
import sys
import tempfile

from tensorflow.keras.models import save_model

import medaka.keras_ext
import medaka.models


def main():
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
    parser = argparse.ArgumentParser(
        description='Convert hdf5 model file to tensorflow model',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'input_model',
        help='Input model filepath.')
    parser.add_argument(
        '--force', action='store_true',
        help='Allow existing model to be overwritten')
    args = parser.parse_args()
    model_name = os.path.splitext(args.input_model)[0]
    out_name = model_name + '.tar.gz'

    if os.path.exists(out_name):
        if not args.force:
            print('{} exists. Use --force to overwrite.'.format(out_name))
            sys.exit(1)
        else:
            os.remove(out_name)

    # get what we need from the old model
    with medaka.datastore.ModelStore(args.input_model) as ms:
        model = ms.load_model(time_steps=None)
        _medaka_meta = {
            'model_function': ms.get_meta('model_function'),
            'label_scheme': ms.get_meta("label_scheme"),
            'feature_encoder': ms.get_meta('feature_encoder')}

    # save new model
    with tempfile.TemporaryDirectory() as tmpdir:
        class Mock:
            epoch_fp = tmpdir
            medaka_meta = _medaka_meta
        mock = Mock()
        save_model(model, tmpdir)
        medaka.keras_ext.ModelMetaCheckpointTF.pack_meta(mock, clean=False)
        shutil.move(tmpdir + '.tar.gz', out_name)
        sys.stdout.write("Model written to {}\n".format(out_name))


if __name__ == '__main__':
    main()

