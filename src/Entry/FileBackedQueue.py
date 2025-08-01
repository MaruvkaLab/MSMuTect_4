import sys, os, time, random


def get_unique_filename():
    return (f".tmp_{os.getpid()}_{time.time()}_{random.randint(1, os.getpid()*4)}") # while not provably random, this is good enough for now


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
