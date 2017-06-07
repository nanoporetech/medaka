

def trim_to_step_multiple(arr, batch_size, window_size):
    """Trim array to be a multiple of batch size
    (keras generator methods require exact numbers of batches)

    :param arr: numpy array
    :param batch_size: int batch size
    :param window_size: int window size
    :returns: array truncated to exact batch multiple
    """
    total_samples = arr.shape[0] - window_size + 1
    remainder = total_samples % batch_size
    arr = arr[:arr.shape[0] - remainder]
    n_batches = (arr.shape[0] - window_size + 1) // batch_size
    return arr, n_batches
