#!/usr/bin/env python3
"""
sqaParallel.py - Universal multiprocessing wrapper for SQA matrix element calculations

This module provides automatic parallelization for computing matrix elements
using the Second Quantization Algebra (SQA) library.

Features:
- Automatic CPU detection
- Hermitian symmetry exploitation
- Optional Slater-Condon selection rules
- Progress reporting with PIDs
- Error handling and verification

Usage example:

```python
import sqaParallel

# Define your compute function (must be at module level for pickling)
def compute_my_element(bra_name, ket_name):
    # Your SQA computation here
    result = ...  # list of SQA terms
    return result

# Define configurations
configs = ['config1', 'config2', 'config3', 'config4']

# Compute matrix with all optimizations
matrix = sqaParallel.compute_matrix(
    compute_func=compute_my_element,
    config_names=configs,
    use_symmetry=True,
    selection_rule_func=None,  # or provide custom function
    n_workers=None  # auto-detect
)

# Access results
result = matrix[('config1', 'config2')]
```

Author: Claude Code
Date: 2025-10-29
"""

import time
import os
from multiprocessing import Pool, cpu_count
from functools import partial


# =========================================================================
# Core Functions
# =========================================================================

def compute_matrix(compute_func,
                  config_names,
                  use_symmetry=True,
                  selection_rule_func=None,
                  n_workers=None,
                  verbose=True):
    """
    Compute full matrix of elements in parallel.

    Parameters
    ----------
    compute_func : callable
        Function with signature: result = compute_func(bra_name, ket_name)
        Returns list of SQA terms.
        MUST be defined at module level for multiprocessing.

    config_names : list of str
        List of configuration names (e.g., ['aabb', 'aacc', ...])

    use_symmetry : bool, default=True
        Use Hermitian symmetry: <I|H|J> = <J|H|I>
        If True, compute only upper triangle + diagonal.

    selection_rule_func : callable or None, default=None
        Optional function to skip zero elements by selection rules.
        Signature: is_zero = selection_rule_func(bra_name, ket_name)
        Returns True if element should be skipped (known to be zero).

    n_workers : int or None, default=None
        Number of parallel workers.
        If None, auto-detect using cpu_count().
        If 0, run sequentially (for debugging).

    verbose : bool, default=True
        Print progress messages.

    Returns
    -------
    matrix : dict
        Dictionary with (bra_name, ket_name) keys and result lists as values.

    statistics : dict
        Dictionary with timing and statistics:
        - 'total_time': Total computation time in seconds
        - 'n_computed': Number of elements computed
        - 'n_skipped': Number of elements skipped by selection rules
        - 'n_copied': Number of elements copied by symmetry
        - 'n_workers': Number of workers used

    Examples
    --------
    >>> def my_compute(bra, ket):
    ...     # Your SQA computation
    ...     return result_terms
    >>>
    >>> matrix, stats = compute_matrix(
    ...     compute_func=my_compute,
    ...     config_names=['aabb', 'aacc', 'bbdd', 'ccdd']
    ... )
    >>>
    >>> print(f"Computed in {stats['total_time']:.1f}s")
    >>> print(f"Speedup: {stats['speedup']:.2f}×")
    """

    # Auto-detect workers
    if n_workers is None:
        n_workers = cpu_count()

    n_configs = len(config_names)
    total_elements = n_configs ** 2

    if verbose:
        print("="*70)
        print("sqaParallel: Matrix Element Computation")
        print("="*70)
        print()
        print(f"System info:")
        print(f"  CPU cores available: {cpu_count()}")
        print(f"  Workers to use: {n_workers if n_workers > 0 else 'sequential'}")
        print(f"  Configurations: {n_configs}")
        print(f"  Total elements: {total_elements}")
        print()

    # Build task list
    tasks = []
    task_id = 0

    if use_symmetry:
        # Upper triangle + diagonal only
        for i, bra_name in enumerate(config_names):
            for j, ket_name in enumerate(config_names):
                if i <= j:
                    tasks.append((bra_name, ket_name, task_id))
                    task_id += 1
    else:
        # All elements
        for bra_name in config_names:
            for ket_name in config_names:
                tasks.append((bra_name, ket_name, task_id))
                task_id += 1

    if verbose:
        print(f"Optimization strategy:")
        if use_symmetry:
            print(f"  Symmetry: ON (compute {len(tasks)} elements, copy {total_elements - len(tasks)})")
        else:
            print(f"  Symmetry: OFF (compute all {len(tasks)} elements)")

        if selection_rule_func:
            # Count how many will be skipped
            n_skippable = sum(1 for bra, ket, _ in tasks
                            if selection_rule_func(bra, ket))
            if n_skippable > 0:
                print(f"  Selection rules: ON (skip {n_skippable} elements)")
            else:
                print(f"  Selection rules: provided but none applicable")
        else:
            print(f"  Selection rules: OFF")
        print()

    # Compute in parallel or sequentially
    start_time = time.time()

    if n_workers == 0:
        # Sequential execution (for debugging)
        if verbose:
            print("Running SEQUENTIALLY (n_workers=0)")
            print()

        results = []
        for task in tasks:
            result = _compute_element_wrapper(task, compute_func, selection_rule_func, verbose)
            results.append(result)
    else:
        # Parallel execution
        if verbose:
            print(f"Running in PARALLEL with {n_workers} workers")
            print("="*70)
            print()

        # Limit workers to number of tasks
        n_workers = min(n_workers, len(tasks))

        # Create wrapper using functools.partial (avoids pickling issues with local functions)
        wrapper = partial(_compute_element_wrapper,
                         compute_func=compute_func,
                         selection_rule_func=selection_rule_func,
                         verbose=verbose)

        with Pool(processes=n_workers) as pool:
            results = pool.map(wrapper, tasks)

    computation_time = time.time() - start_time

    # Build full matrix
    matrix = {}
    n_computed = 0
    n_skipped = 0

    # Add computed/skipped elements
    for (bra_name, ket_name, result, elapsed, skipped) in results:
        matrix[(bra_name, ket_name)] = result
        if skipped:
            n_skipped += 1
        else:
            n_computed += 1

    # Add symmetry copies
    n_copied = 0
    if use_symmetry:
        for i, bra_name in enumerate(config_names):
            for j, ket_name in enumerate(config_names):
                if i > j:  # Lower triangle
                    matrix[(bra_name, ket_name)] = matrix[(ket_name, bra_name)]
                    n_copied += 1

    total_time = time.time() - start_time

    # Statistics
    statistics = {
        'total_time': total_time,
        'computation_time': computation_time,
        'n_computed': n_computed,
        'n_skipped': n_skipped,
        'n_copied': n_copied,
        'n_workers': n_workers,
        'total_elements': total_elements
    }

    # Calculate speedup estimate
    if n_skipped > 0 or n_copied > 0:
        # Sequential time = time for all unique elements without optimizations
        estimated_sequential = computation_time * total_elements / (n_computed + n_skipped)
        statistics['speedup'] = estimated_sequential / computation_time
    else:
        # Just parallel speedup (no selection rules or symmetry)
        statistics['speedup'] = n_workers

    if verbose:
        print()
        print("="*70)
        print("COMPUTATION STATISTICS")
        print("="*70)
        print(f"Workers:            {n_workers}")
        print(f"Elements computed:  {n_computed}")
        if n_skipped > 0:
            print(f"Elements skipped:   {n_skipped} (selection rules)")
        if n_copied > 0:
            print(f"Elements copied:    {n_copied} (symmetry)")
        print(f"Total elements:     {total_elements}")
        print(f"Computation time:   {computation_time:.1f}s ({computation_time/60:.2f} min)")
        print(f"Total time:         {total_time:.1f}s ({total_time/60:.2f} min)")
        print(f"Speedup estimate:   {statistics['speedup']:.2f}×")
        print()

    return matrix, statistics


def _compute_element_wrapper(task, compute_func, selection_rule_func, verbose):
    """
    Internal wrapper for parallel execution.

    Parameters
    ----------
    task : tuple
        (bra_name, ket_name, task_id)
    compute_func : callable
        Function to compute matrix element
    selection_rule_func : callable or None
        Function to check if element is zero by selection rules
    verbose : bool
        Print progress messages

    Returns
    -------
    result : tuple
        (bra_name, ket_name, result_terms, elapsed_time, skipped)
    """
    bra_name, ket_name, task_id = task
    pid = os.getpid()

    # Check selection rules first
    if selection_rule_func and selection_rule_func(bra_name, ket_name):
        if verbose:
            print(f"[Task {task_id}] PID {pid}: <{bra_name}|H|{ket_name}> SKIPPED (selection rule)", flush=True)
        return (bra_name, ket_name, [], 0.0, True)

    # Compute normally
    if verbose:
        print(f"[Task {task_id}] PID {pid}: Starting <{bra_name}|H|{ket_name}>", flush=True)

    start_time = time.time()

    try:
        result = compute_func(bra_name, ket_name)
        elapsed = time.time() - start_time

        if verbose:
            if len(result) == 0:
                print(f"[Task {task_id}] PID {pid}: <{bra_name}|H|{ket_name}>: = 0 [{elapsed:.1f}s]", flush=True)
            else:
                print(f"[Task {task_id}] PID {pid}: <{bra_name}|H|{ket_name}>: {len(result)} terms [{elapsed:.1f}s]", flush=True)

        return (bra_name, ket_name, result, elapsed, False)

    except Exception as e:
        elapsed = time.time() - start_time
        if verbose:
            print(f"[Task {task_id}] PID {pid}: <{bra_name}|H|{ket_name}>: ERROR: {e} [{elapsed:.1f}s]", flush=True)
        return (bra_name, ket_name, None, elapsed, False)


# =========================================================================
# Selection Rule Helpers
# =========================================================================

def slater_condon_diff(config_a, config_b):
    """
    Compute diff = number of spin-orbitals that differ between two configurations.

    This is the standard definition for Slater-Condon rules:
        diff = |A \ B| = |B \ A|
    where A and B are SETS of occupied spin-orbitals.

    Parameters
    ----------
    config_a : list or set
        List/set of orbital names in configuration A
    config_b : list or set
        List/set of orbital names in configuration B

    Returns
    -------
    diff : int
        Number of spin-orbitals that differ

    Examples
    --------
    >>> slater_condon_diff(
    ...     ['a_alpha', 'a_beta', 'b_alpha', 'b_beta'],
    ...     ['a_alpha', 'a_beta', 'c_alpha', 'c_beta']
    ... )
    2
    """
    set_a = set(config_a) if not isinstance(config_a, set) else config_a
    set_b = set(config_b) if not isinstance(config_b, set) else config_b

    # diff = number of orbitals in A but not in B
    # (equals number in B but not in A)
    diff = len(set_a - set_b)

    return diff


def make_selection_rule_two_electron(configurations, threshold=3):
    """
    Create selection rule function for two-electron operators.

    For two-electron operators, Slater-Condon rules give:
        <I|H_{2e}|J> = 0  if  diff(I,J) >= threshold

    Standard threshold is 3 for two-electron operators.

    Parameters
    ----------
    configurations : dict
        Dictionary mapping config names to orbital lists
        Example: {'aabb': ['a_alpha', 'a_beta', 'b_alpha', 'b_beta'], ...}

    threshold : int, default=3
        Diff threshold above which elements are zero

    Returns
    -------
    selection_rule_func : callable
        Function with signature: is_zero = func(bra_name, ket_name)

    Examples
    --------
    >>> configs = {
    ...     'aabb': ['a_alpha', 'a_beta', 'b_alpha', 'b_beta'],
    ...     'ccdd': ['c_alpha', 'c_beta', 'd_alpha', 'd_beta']
    ... }
    >>> rule = make_selection_rule_two_electron(configs)
    >>> rule('aabb', 'ccdd')  # diff=4, should be zero
    True
    >>> rule('aabb', 'aabb')  # diff=0, non-zero
    False
    """
    def is_zero_by_selection_rule(bra_name, ket_name):
        bra_orbs = configurations[bra_name]
        ket_orbs = configurations[ket_name]
        diff = slater_condon_diff(bra_orbs, ket_orbs)
        return diff >= threshold

    return is_zero_by_selection_rule


# =========================================================================
# Verification Helpers
# =========================================================================

def verify_symmetry(matrix, config_names, verbose=True):
    """
    Verify Hermitian symmetry: <I|H|J> = <J|H|I>

    Parameters
    ----------
    matrix : dict
        Matrix dictionary from compute_matrix()
    config_names : list
        List of configuration names
    verbose : bool, default=True
        Print verification results

    Returns
    -------
    all_symmetric : bool
        True if all symmetry checks pass

    Examples
    --------
    >>> matrix, stats = compute_matrix(...)
    >>> verify_symmetry(matrix, config_names)
    ✓ All symmetry checks passed!
    True
    """
    if verbose:
        print("="*70)
        print("Symmetry Verification")
        print("="*70)
        print()

    all_symmetric = True

    for i, bra in enumerate(config_names):
        for j, ket in enumerate(config_names):
            if i < j:  # Upper triangle only
                H_IJ = matrix[(bra, ket)]
                H_JI = matrix[(ket, bra)]

                # Check if same
                if len(H_IJ) == len(H_JI):
                    # Check each term
                    terms_match = True
                    for term_ij, term_ji in zip(H_IJ, H_JI):
                        if str(term_ij) != str(term_ji):
                            terms_match = False
                            break

                    if terms_match:
                        if verbose:
                            print(f"✓ <{bra}|H|{ket}> = <{ket}|H|{bra}>")
                    else:
                        if verbose:
                            print(f"✗ <{bra}|H|{ket}> ≠ <{ket}|H|{bra}> (terms differ)")
                        all_symmetric = False
                else:
                    if verbose:
                        print(f"✗ <{bra}|H|{ket}> ≠ <{ket}|H|{bra}> (length: {len(H_IJ)} vs {len(H_JI)})")
                    all_symmetric = False

    if verbose:
        print()
        if all_symmetric:
            print("✓ All symmetry checks passed!")
        else:
            print("✗ Symmetry violations detected!")
        print()

    return all_symmetric


# =========================================================================
# Example Usage
# =========================================================================

if __name__ == '__main__':
    print(__doc__)
    print()
    print("This is a library module. Import it in your script:")
    print()
    print("  import sqaParallel")
    print()
    print("See docstrings for usage examples.")
