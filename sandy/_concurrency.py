"""
Internal concurrency helpers for spawn-safe, cross-platform pools.
"""



"""
`_G_ENDF6`, `_G_PENDF`, `_G_RDD`, `_G_NFPY` “globals” are process‑local, not shared memory.
Each worker process gets its own copy.
They’re initialized once per worker by the pool’s initializer=... (not per task),
and then reused for all tasks executed by that worker until the worker exits.

The goal of the init_*_cache functions is to load the large, read‑only “nominal” 
objects (ENDF‑6/PENDF, RDD, NFY) once per worker process.
"""
# Per-process caches (filled via initializer)
_G_ENDF6 = None      # dict (nominal ENDF6)
_G_PENDF = None      # dict (nominal PENDF)
_G_RDD = None        # DecayData.data (dict-like)
_G_NFPY = None       # Fy.data (DataFrame)



def spawn_ctx():
    """Always use a spawn context (Windows default; explicit elsewhere)."""
    from multiprocessing import get_context
    return get_context("spawn")



def init_xs_cache(endf6_dict, pendf_dict):
    """
    Initialize per‑process, read‑only caches for ENDF6 and PENDF data
    used by XS perturbation workers.

    This function is designed to be passed as the `initializer` argument of a
    `concurrent.futures.ProcessPoolExecutor` (or `multiprocessing.Pool`).
    The executor calls it **once per worker process**, right after the worker
    starts—so large inputs are set into per‑process globals only once and then
    reused by all tasks executed in that worker.

    Why this matters
    ----------------
    - **Windows/macOS (spawn)**: worker processes start from a fresh
      interpreter. Without a per‑process cache, each task would repeatedly
      pickle and ship large data structures to every worker call (slow and
      memory‑intensive). Initializing a cache **once per worker** avoids that
      overhead.
    - **Linux (often fork)**: even though forking can be copy‑on‑write
      efficient, using an initializer keeps behavior consistent (and is still
      beneficial if you explicitly choose `spawn`).

    Parameters
    ----------
    endf6_dict : dict
        The nominal ENDF‑6 content for the evaluation being perturbed.
        This must be the plain `dict` you would normally find in
        `Endf6(...).data` (i.e., *not* an `Endf6` instance), so it can be
        pickled once and cached in the worker.
    pendf_dict : dict
        The nominal PENDF content (again, plain `dict`, not an `Endf6`
        object). It serves as the starting point for cross‑section
        perturbations before any Doppler broadening / reconstruction steps
        performed per sample.

    Side Effects
    ------------
    Sets (or replaces) two module‑level globals in the worker process:
    - `_G_ENDF6`: the cached ENDF‑6 dict
    - `_G_PENDF`: the cached PENDF dict

    These globals are **process‑local**: every worker process has its own
    copy. They are intentionally **read‑only** for tasks (do not mutate them).

    """
    # Implementation (kept minimal on purpose):
    global _G_ENDF6, _G_PENDF
    _G_ENDF6 = endf6_dict
    _G_PENDF = pendf_dict



def init_rdd_cache(endf6_dict, rdd_data):
    """
    Initializer run ONCE per worker process. Caches nominal ENDF6 and RDD dicts
    into per-process globals. Mirrors the XS initializer pattern.
    """
    global _G_ENDF6, _G_RDD
    _G_ENDF6 = endf6_dict
    _G_RDD = rdd_data



def init_fy_cache(endf6_dict, nfpy_data):
    global _G_ENDF6, _G_NFPY
    _G_ENDF6 = endf6_dict
    _G_NFPY = nfpy_data



# Thin task wrappers that use the caches. We import lazily inside functions
# to avoid circular imports and keep import time low.
def task_xs(
        ismp,
        *,
        pxs=None,
        pnu=None,
        pchi=None,
        verbose: bool = False,
        to_ace: bool =False,
        to_file: bool = False,
        ace_kws=None,
        ):
    """
    Execute a single **XS/nubar/chi perturbation** task using per‑process caches.

    This wrapper is the **public concurrency entry point** for ENDF‑6 neutron
    data perturbations (XS, nubar, chi) when running inside a spawned worker
    process. It forwards a **small, explicit set of arguments** to the internal
    domain worker (`_endf6_perturb_worker`) while obtaining the large, nominal
    inputs (ENDF‑6/PENDF) from **per‑process globals** that were previously
    initialized once by `init_xs_cache`.

    Why this wrapper exists
    -----------------------
    1) **Performance on Windows/macOS (spawn)**:
       - Worker processes start from a fresh interpreter and everything passed to
         them is pickled. If we were to ship the large ENDF‑6/PENDF structures
         with *every* task, the overhead would be severe.
       - `task_xs` relies on per‑process caches (`_G_ENDF6`, `_G_PENDF`)
         populated by `init_xs_cache` **once per worker**; tasks then only pass
         the *small* per‑sample payloads (e.g., a single perturbation dataframe).

    2) **Safety & clarity**:
       - The **explicit signature** (`pxs`, `pnu`, `pchi`, …) documents exactly
         what the task accepts and keeps a stable boundary between the
         concurrency layer and the domain worker.
       - It also prevents accidental leakage of unrelated `**kwargs`—we forward
         only what we intend to support.

    3) **Process‑local isolation**:
       - The caches are **process‑local**, not shared memory. Each worker has
         its own copy of the nominal ENDF‑6/PENDF. Treat these as **read‑only**
         inside tasks; the domain worker already copies before mutating.

    Parameters
    ----------
    Same as in :func:`sandy._workers._endf6_perturb_worker`.

    Returns
    -------
    Same as in :func:`sandy._workers._endf6_perturb_worker`.

    Preconditions
    -------------
    - The worker process **must** have run `init_xs_cache(endf6_dict, pendf_dict)`
      as its `initializer` when the pool/executor started; otherwise the caches
      `_G_ENDF6`/`_G_PENDF` will be `None` and the call will fail.

    Why not call the internal worker directly?
    ------------------------------------------
    - :func:`sandy._workers._endf6_perturb_worker` expects **nominal ENDF‑6/PENDF**
      each time. Passing those from the parent for every task would cause **heavy pickling**.
    - `task_xs` enforces the **small‑payload pattern** and ensures that the
      nominal inputs come from the **per‑process cache** instead.

    See also
    --------
    :func:`sandy._concurrency.init_xs_cache`     : initializer that sets `_G_ENDF6` and `_G_PENDF` once per worker
    :func:`sandy._concurrency.task_rdd`          : analogous wrapper for radioactive decay perturbations
    :func:`sandy._concurrency.task_fy`           : analogous wrapper for fission yield perturbations
    :func:`sandy._workers._endf6_perturb_worker` : internal domain worker that performs the actual changes
    """
    # --- implementation: thin wrapper over the domain worker, using per-process caches
    from ._workers import _endf6_perturb_worker  # local import to avoid cycles

    return _endf6_perturb_worker(
        endf6=_G_ENDF6,
        pendf=_G_PENDF,
        ismp=ismp,
        pxs=pxs,
        pnu=pnu,
        pchi=pchi,
        verbose=verbose,
        to_ace=to_ace,
        to_file=to_file,
        ace_kws=ace_kws,
    )



def task_rdd(
        ismp,
        phl,
        pde,
        pbr,
        *,
        verbose: bool = False,
        to_file: bool = False,
        ):
    """
    Per-sample task. Uses cached nominal dicts and per-sample perturbation columns.
    """
    from ._workers import _rdd_perturb_worker

    # Delegate to existing low-level worker
    return _rdd_perturb_worker(
        _G_ENDF6,
        _G_RDD,
        phl,
        pde,
        pbr,
        ismp,
        verbose=verbose,
        to_file=to_file,
    )



def task_fy(
        ismp,
        pfy,
        *,
        verbose: bool = False,
        to_file: bool = False,
        ):
    """
    Per-sample task. Uses cached nominal dicts and per-sample perturbation columns.
    """
    from ._workers import _fy_perturb_worker
    
    # Delegate to existing low-level worker
    return _fy_perturb_worker(
        endf6=_G_ENDF6,
        fy=_G_NFPY,
        pfy=pfy,
        ismp=ismp,
        verbose=verbose,
        to_file=to_file,
        )