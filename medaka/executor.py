"""Wrappers around standard library multiprocessing executors."""
import concurrent.futures
import multiprocessing
import threading


class _Executor(object):

    def submit(self, fn, *args, **kwargs):
        self.semaphore.acquire()
        future = super().submit(fn, *args, **kwargs)
        future.add_done_callback(self._release)

        return future

    def _release(self, future):
        self.semaphore.release()


class ProcessPoolExecutor(_Executor, concurrent.futures.ProcessPoolExecutor):
    """Extends `ProcessPoolExecutor` by limiting simultaneous work items."""

    def __init__(self, max_items, **kwargs):
        """Initialize a process pool.

        :param max_items: maximum number of simultaneous work items.
            Calls to `.submit` will block if there are too many unprocessed
            items.
        :param kwargs: key-word arguments to `ProcessPoolExecutor`.
        """
        super().__init__(**kwargs)
        self.semaphore = multiprocessing.BoundedSemaphore(max_items)


class ThreadPoolExecutor(_Executor, concurrent.futures.ThreadPoolExecutor):
    """Extends `ThreadPoolExecutor` by limiting simultaneous work items."""

    def __init__(self, max_items, **kwargs):
        """Initialize a thread pool.

        :param max_items: maximum number of simultaneous work items.
            Calls to `.submit` will block if there are too many unprocessed
            items.
        :param kwargs: key-word arguments to `ThreadPoolExecutor`.
        """
        super().__init__(**kwargs)
        self.semaphore = threading.BoundedSemaphore(max_items)
