import sys, os, time, random


TMP_FILE_PREFIX = ".msmutect_tmp_file_"
RUN_ID_ENV_VAR = "MSMUTECT_RUN_ID"


# Process-local random run id, used when the environment doesn't provide one
# (standalone use of FileBackedQueue). Generated once so every temp file from this
# process shares it and can be cleaned up together, and random rather than pid-based
# so it stays unique even across pid reuse.
_FALLBACK_RUN_ID = f"{os.getpid()}_{random.randint(0, 2 ** 31)}"


def get_run_id() -> str:
    # Stable per-run token set by the parent process before any workers fork
    # (see main.run_msmutect). Forked/spawned workers inherit it through the
    # environment, so every temp file produced by one run shares the same run id
    # and can be cleaned up precisely. Falls back to a process-local random id
    # when unset.
    return os.environ.get(RUN_ID_ENV_VAR) or _FALLBACK_RUN_ID


def get_unique_filename():
    # run id (shared, for scoped cleanup) + timestamp + a random suffix. The random
    # component makes the name unique independent of the pid, so concurrent runs
    # sharing an output directory never collide.
    return f"{TMP_FILE_PREFIX}{get_run_id()}_{time.time()}_{random.randint(0, 2 ** 31)}"


class FileBackedQueue:
    """
    Will write items in itself to a file when it passes a certain size
    """
    def __init__(self, out_file_dir: str = "", max_memory: int = 10**7):
        out_file = get_unique_filename()
        self.out_file_path = os.path.join(out_file_dir, out_file)
        self.out_file = open(f"{self.out_file_path}", 'w+')
        self.max_memory = max_memory
        self.queue = []
        self.queue_size = 0

    def append(self, item: str):
        self.queue.append(item)
        self.queue_size+=sys.getsizeof(item)
        if self.queue_size > self.max_memory:
            self.flush()

    def flush(self):
        self.out_file.write("\n".join(self.queue))
        self.out_file.write("\n")
        self.queue = []
        self.queue_size = 0

    def close(self):
        # should be called when all results have been written
        self.flush()
        self.out_file.close()
        del self.out_file

    def delete_backing_file(self):
        os.remove(self.out_file_path)

    def __del__(self):
        try:
            self.out_file.close()
        except AttributeError: # already deleted
            pass
