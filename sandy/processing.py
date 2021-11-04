import multiprocessing as mp

__all__ = [
    ]
__author__ = "Luca Fiorito"


def run_parallel(foo, entries, processes, *args, **kwargs):
    pool = mp.Pool(processes=processes)
    outs = {x: pool.apply_async(foo, (x, *args), kwargs) for x in entries}
    outs = {x: out.get() for x, out in outs.items()}
    pool.close()
    pool.join()
    return outs
