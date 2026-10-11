import multiprocessing.pool, time, os, shutil
from typing import List
from collections import namedtuple
from multiprocessing import Pool

from src.Entry.FileBackedQueue import FileBackedQueue
from src.GenomicUtils.LocusFile import LociManager

Chunk = namedtuple("Chunk", ["start", "end"])

def minimum_chunk_size() -> int:
    return 20_000

def get_noise_table_path() -> str:
    script_dir = os.path.dirname(os.path.realpath(__file__))
    noise_table_path = script_dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'data/noise_table.csv'
    return os.path.abspath(noise_table_path)


def get_batch_sizes(total_batch_size: int, regular_size: int) -> List[int]:
    batch_sizes = [regular_size for _ in range(total_batch_size // regular_size)]
    remainder = total_batch_size % regular_size
    if remainder != 0:
        batch_sizes.append(remainder)
    return batch_sizes


def extract_results(results: List[multiprocessing.pool.ApplyResult]) -> List[FileBackedQueue]:
    # extracts results from multiproccessing
    combined = [result.get() for result in results]
    return combined


def write_results(output_prefix: str, results: List[str], header):
    with open(f"{output_prefix}.tsv", 'w+') as output_file:
        output_file.write(header)
        output_file.write("\n")
        output_file.write("\n".join(results))


def write_queues_results(output_prefix: str, results: List[FileBackedQueue], header: str):
    with open(f"{output_prefix}.tsv", 'w+') as output_file:
        output_file.write(header)
        output_file.write("\n")
        for r in results:
            with open(r.out_file_path, 'r') as q_file:
                shutil.copyfileobj(q_file, output_file)
                # output_file.write("\n")
            r.delete_backing_file()


def run_single_threaded(batch_function, args: list, loci_iterator: LociManager, total_batch_size: int, result_dir: str) -> list:
    """
    runs batch fuction without invoking pool to save performance (serialization, etc.)
    """
    results = []
    batch_sizes = get_batch_sizes(total_batch_size, minimum_chunk_size())
    for batch in batch_sizes:
        current_loci = loci_iterator.get_batch(batch)
        results.append(batch_function(*([current_loci] + args + [result_dir])))
    return results


def run_batch(batch_function, args: list, loci_iterator: LociManager, total_batch_size: int, cores: int, result_dir: str) -> List[FileBackedQueue]:
    """
    :param batch_function: function to run on given loci. First argument must be list of loci
    :param args: other args to feed function
    :return: results from given function
    """
    results = []
    if cores == 1:
        return run_single_threaded(batch_function, args, loci_iterator, total_batch_size, result_dir)
    with Pool(processes=cores) as threads:
        batch_sizes = get_batch_sizes(total_batch_size, minimum_chunk_size())

        for batch in batch_sizes:
            num_active_processes = sum([1 for p in results if not p.ready()]) # how many processes are actually running
            while num_active_processes == cores:
                time.sleep(1)  # long wait time to avoid wasting processing power
                num_active_processes = sum([1 for p in results if not p.ready()])
            current_loci = loci_iterator.get_batch(batch)
            results.append(threads.apply_async(batch_function,
                                           args=([current_loci] + args + [result_dir])))


        # for i, batch in enumerate(batch_sizes):
        #     current_loci = loci_iterator.get_batch(batch)
        #     results.append(threads.apply_async(batch_function,
        #                                        args=([current_loci] + args + [result_dir])))
        threads.close()
        threads.join()
    return extract_results(results)
