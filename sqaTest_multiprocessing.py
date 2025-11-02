#!/usr/bin/env python3
"""
sqaTest_multiprocessing.py - Multiprocessing Performance Test

Tests the sqaParallel module for computing matrix elements in parallel.

This test demonstrates ~4× speedup on 4-core systems by computing
4 diagonal matrix elements simultaneously.

Expected runtime:
- Sequential: ~2 minutes
- Parallel (4 cores): ~30 seconds
- Speedup: ~4×

Usage:
    python3 sqaTest_multiprocessing.py
"""
import sys
import time
sys.path.insert(0, '.')

import secondQuantizationAlgebra as sqa

# Disable verbose output
sqa.options.verbose = False

print("="*70)
print("SQA Multiprocessing Test")
print("="*70)
print()

# =========================================================================
# Setup: Define 4-electron system configurations
# =========================================================================
a_alpha = sqa.index('a_alpha', [], False)
a_beta = sqa.index('a_beta', [], False)
b_alpha = sqa.index('b_alpha', [], False)
b_beta = sqa.index('b_beta', [], False)
c_alpha = sqa.index('c_alpha', [], False)
c_beta = sqa.index('c_beta', [], False)
d_alpha = sqa.index('d_alpha', [], False)
d_beta = sqa.index('d_beta', [], False)

all_orbs_dict = {
    'a_alpha': a_alpha, 'a_beta': a_beta,
    'b_alpha': b_alpha, 'b_beta': b_beta,
    'c_alpha': c_alpha, 'c_beta': c_beta,
    'd_alpha': d_alpha, 'd_beta': d_beta
}

configurations = {
    'aabb': ['a_alpha', 'a_beta', 'b_alpha', 'b_beta'],
    'aacc': ['a_alpha', 'a_beta', 'c_alpha', 'c_beta'],
    'bbdd': ['b_alpha', 'b_beta', 'd_alpha', 'd_beta'],
    'ccdd': ['c_alpha', 'c_beta', 'd_alpha', 'd_beta']
}

# =========================================================================
# Helper functions
# =========================================================================
def build_ket(orb_names):
    ops = [sqa.creOp(all_orbs_dict[name]) for name in orb_names]
    return sqa.term(1.0, [], ops)

def build_bra(orb_names):
    ops = [sqa.desOp(all_orbs_dict[name]) for name in reversed(orb_names)]
    return sqa.term(1.0, [], ops)

def same_spin(orb1, orb2):
    both_alpha = orb1.name.endswith('alpha') and orb2.name.endswith('alpha')
    both_beta = orb1.name.endswith('beta') and orb2.name.endswith('beta')
    return both_alpha or both_beta

def build_H2e_for_orbitals(orb_names):
    orbs = [all_orbs_dict[name] for name in orb_names]
    ops_cre = [sqa.creOp(orb) for orb in orbs]
    ops_des = [sqa.desOp(orb) for orb in orbs]

    H2e_terms = []
    n = len(orbs)

    for i in range(n):
        for j in range(i+1, n):
            for k in range(n):
                for l in range(n):
                    if k == l:
                        continue
                    is_exchange = (k == j and l == i)
                    if is_exchange and not same_spin(orbs[i], orbs[j]):
                        continue

                    g2 = sqa.tensor('g2', [orbs[i], orbs[j], orbs[k], orbs[l]], [])
                    term = sqa.term(0.5, [], [g2, ops_cre[i], ops_cre[j],
                                               ops_des[l], ops_des[k]])
                    H2e_terms.append(term)

    return H2e_terms

def compute_diagonal_element(config_name):
    """Compute <config|H_{2e}|config> matrix element."""
    orb_names = configurations[config_name]

    H2e_terms = build_H2e_for_orbitals(orb_names)
    bra = build_bra(orb_names)
    ket = build_ket(orb_names)

    result_terms = []

    for H_term in H2e_terms:
        product1 = sqa.multiplyTerms(H_term, ket)
        h_ket = sqa.normalOrder(product1)

        for h_ket_term in h_ket:
            product2 = sqa.multiplyTerms(bra, h_ket_term)
            result = sqa.normalOrder(product2)

            for term in result:
                has_ops = any(isinstance(t, (sqa.creOp, sqa.desOp)) for t in term.tensors)
                if not has_ops:
                    single_term_list = [term]
                    sqa.contractDeltaFuncs_complete(single_term_list, zero_remaining=True)

                    if len(single_term_list) > 0:
                        for contracted_term in single_term_list:
                            has_kdelta = any(t.name == 'kdelta' for t in contracted_term.tensors)
                            if not has_kdelta and abs(contracted_term.numConstant) > 1e-10:
                                result_terms.append(contracted_term)

    if len(result_terms) > 0:
        sqa.combineTerms(result_terms)

    return result_terms

# =========================================================================
# Test 1: Sequential Execution
# =========================================================================
print("Test 1: Sequential Execution")
print("-" * 70)

config_names = ['aabb', 'aacc', 'bbdd', 'ccdd']

seq_start = time.time()
seq_results = {}

for config_name in config_names:
    result = compute_diagonal_element(config_name)
    seq_results[config_name] = result
    print(f"  <{config_name}|H|{config_name}>: {len(result)} terms")

seq_time = time.time() - seq_start

print(f"\nSequential time: {seq_time:.1f}s ({seq_time/60:.2f} min)")
print()

# =========================================================================
# Test 2: Parallel Execution
# =========================================================================
print("Test 2: Parallel Execution")
print("-" * 70)

try:
    import sqaParallel

    # Compute in parallel
    par_start = time.time()

    matrix, stats = sqaParallel.compute_matrix(
        compute_func=compute_diagonal_element,
        config_names=config_names,
        use_symmetry=False,  # Only diagonal elements
        selection_rule_func=None,
        n_workers=None,  # auto-detect
        verbose=False
    )

    par_time = time.time() - par_start

    print(f"Parallel time: {par_time:.1f}s ({par_time/60:.2f} min)")
    print(f"Workers used: {stats['n_workers']}")
    print()

    # =========================================================================
    # Results
    # =========================================================================
    print("="*70)
    print("Results")
    print("="*70)
    print(f"Sequential time: {seq_time:.1f}s ({seq_time/60:.2f} min)")
    print(f"Parallel time:   {par_time:.1f}s ({par_time/60:.2f} min)")
    print(f"Speedup:         {seq_time/par_time:.2f}×")
    print()

    # Verify results match
    print("Verification:")
    all_match = True
    for config_name in config_names:
        seq = seq_results[config_name]
        par = matrix[(config_name, config_name)]

        if len(seq) == len(par):
            terms_match = all(str(s) == str(p) for s, p in zip(seq, par))
            if terms_match:
                print(f"  ✓ {config_name}: {len(seq)} terms match")
            else:
                print(f"  ✗ {config_name}: Terms differ!")
                all_match = False
        else:
            print(f"  ✗ {config_name}: Different lengths ({len(seq)} vs {len(par)})")
            all_match = False

    print()
    if all_match:
        print("✓ All results match! Multiprocessing works correctly.")
        print()
        print("="*70)
        print("TEST PASSED")
        print("="*70)
    else:
        print("✗ Results don't match! Check implementation.")
        print()
        print("="*70)
        print("TEST FAILED")
        print("="*70)
        sys.exit(1)

except ImportError:
    print("WARNING: sqaParallel module not found")
    print("Parallel test skipped")
    print()
    print("="*70)
    print("TEST SKIPPED (sqaParallel not available)")
    print("="*70)
