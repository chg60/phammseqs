"""Wrapper around a subset of the multiprocessing built-in library
to make running multiple processes easy. Load balancing is left to
the caller."""

import multiprocessing

CPUS = multiprocessing.cpu_count()


def show_progress(current, end, width=50):
    """Format an inline-updatable progress bar and print it to the
    console.

    :param current: current step (1 through n)
    :type current: int
    :param end: number of steps (n)
    :type end: int
    :param width: character width for the progressbar
    :type width: int
    """
    width = int(width)

    completion = int(float(current) / end * 100)
    per_symbol = int(100 / width)
    num_filled = int(completion / per_symbol)
    num_padded = width - num_filled

    bar = f"[{'#' * num_filled}{' ' * num_padded}] {completion}%"

    if completion < 100:
        print("\r" + bar, end="")
    else:
        print("\r" + bar)


def parallelize(inputs, cpus, task, verbose=False):
    """Parallelize some task on an input list across the specified
    number of CPU cores.

    :param inputs: arguments to call `task` with
    :type inputs: list
    :param cpus: number of processor cores to use
    :type cpus: int
    :param task: name of the function to run_clustalo
    :type task: function
    :param verbose: updating progress bar output?
    :return: results
    """
    results = []

    # Don't do any work if there are no inputs
    if len(inputs) == 0:
        return results

    # User requested some number of cpus - make sure it's sane
    if cpus < 1 or cpus > CPUS:
        cpus = min([CPUS, len(inputs)])

    tasks = []
    for item in inputs:
        if not isinstance(item, tuple):
            item = (item,)
        tasks.append((task, item))

    # Start working on the jobs
    results = start_processes(tasks, cpus, verbose)

    return results


def worker(input_queue, output_queue):
    """Worker function to run_clustalo jobs from the input queue until a stop
    signal is reached.

    :param input_queue: threadsafe queue to pull jobs from
    :type input_queue: mp.Queue
    :param output_queue: threadsafe queue to add results to
    :type output_queue: mp.Queue
    """
    for func, args in iter(input_queue.get, 'STOP'):
        output_queue.put(func(*args))


def start_processes(inputs, cpus, verbose=False):
    """Spool processes and use them to run_clustalo jobs.

    :param inputs: jobs to run_clustalo
    :param cpus: optimized number of processors
    :param verbose: updating progress bar output?
    :return: results
    """
    job_q = multiprocessing.Queue()
    result_q = multiprocessing.Queue()

    # Calculate how often to place a show_progress() job
    n_tasks = len(inputs)
    frequency = max([1, len(inputs)//100])
    for i, job in enumerate(inputs):
        job_q.put(job)
        if verbose and (i % frequency == 0 or i == len(inputs) - 1):
            job_q.put((show_progress, (i+1, len(inputs))))
            n_tasks += 1

    # Set up workers and put a 'STOP' signal at the end of job_q for each
    worker_pool = list()
    for i in range(cpus):
        worker_pool.append(
            multiprocessing.Process(target=worker, args=(job_q, result_q)))
        job_q.put('STOP')

    # Ready... set... go!
    [w.start() for w in worker_pool]

    # Grab results from result_q
    results = []
    for _ in range(n_tasks):
        result = result_q.get()
        if result is not None:
            results.append(result)

    [w.join() for w in worker_pool]

    return results
