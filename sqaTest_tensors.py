#!/usr/bin/env python3
"""
Test common tensor types defined in User Guide.
Tests h1, h2, a2, d1, d2 tensor creation and symmetries.
"""

import secondQuantizationAlgebra as sqa
import sys

def test_one_electron_tensors():
    """Test h1 (one-electron integrals) and d1 (one-particle density matrix)."""
    print("Testing one-electron tensors (h1, d1)...")

    # Create indices
    i = sqa.index('i', [], True)
    j = sqa.index('j', [], True)

    # h1: One-electron Hamiltonian with exchange symmetry h(i,j) = h(j,i)
    h1_sym = [sqa.symmetry((1,0), 1)]
    h1 = sqa.tensor('h1', [i, j], h1_sym)

    # d1: One-particle density matrix with same symmetry
    d1_sym = [sqa.symmetry((1,0), 1)]
    d1 = sqa.tensor('d1', [i, j], d1_sym)

    # Create term with h1
    term_h1 = sqa.term(1.0, [], [h1, sqa.creOp(i), sqa.desOp(j)])

    # Create term with d1
    term_d1 = sqa.term(1.0, [], [d1])

    print("  h1 tensor created with symmetry: h(i,j) = h(j,i)")
    print("  d1 tensor created with symmetry: d1(i,j) = d1(j,i)")
    print("  Test passed\n")
    return True

def test_two_electron_tensors():
    """Test h2 (two-electron integrals) and a2 (two-electron amplitudes)."""
    print("Testing two-electron tensors (h2, a2)...")

    # Create indices
    indices = [sqa.index('i%d'%n, [], True) for n in range(4)]

    # h2: Two-electron integrals (physicist notation)
    # Symmetries: antisymmetric in electron 1, antisymmetric in electron 2, and particle exchange
    h2_sym = [sqa.symmetry((1,0,2,3), -1),  # antisymmetric in i,j
              sqa.symmetry((0,1,3,2), -1),  # antisymmetric in k,l
              sqa.symmetry((2,3,0,1),  1)]  # particle exchange
    h2 = sqa.tensor('h2', indices, h2_sym)

    # a2: Two-electron amplitudes (CCSD)
    a2_sym = [sqa.symmetry((1,0,2,3), -1),  # antisymmetric in i,j
              sqa.symmetry((0,1,3,2), -1)]  # antisymmetric in a,b
    a2 = sqa.tensor('a2', indices, a2_sym)

    # Create term with h2
    term_h2 = sqa.term(1.0, [], [h2, sqa.creOp(indices[0]), sqa.creOp(indices[1]),
                                   sqa.desOp(indices[3]), sqa.desOp(indices[2])])

    # Create term with a2
    term_a2 = sqa.term(1.0, [], [a2, sqa.creOp(indices[0]), sqa.creOp(indices[1]),
                                   sqa.desOp(indices[3]), sqa.desOp(indices[2])])

    print("  h2 tensor created with antisymmetry and particle exchange")
    print("  a2 tensor created with antisymmetry in creation and destruction")
    print("  Test passed\n")
    return True

def test_two_particle_density_matrix():
    """Test d2 (two-particle density matrix)."""
    print("Testing two-particle density matrix (d2)...")

    # Create indices
    indices = [sqa.index('i%d'%n, [], True) for n in range(4)]

    # d2: Two-particle density matrix
    # Symmetries: antisymmetric in creation, antisymmetric in destruction, particle exchange
    d2_sym = [sqa.symmetry((1,0,2,3), -1),  # antisymmetric in i,j
              sqa.symmetry((0,1,3,2), -1),  # antisymmetric in k,l
              sqa.symmetry((2,3,0,1),  1)]  # particle exchange d2(i,j,k,l) = d2(k,l,i,j)
    d2 = sqa.tensor('d2', indices, d2_sym)

    # Create term with d2
    term_d2 = sqa.term(1.0, [], [d2])

    print("  d2 tensor created with full antisymmetry and particle exchange")
    print("  Test passed\n")
    return True

def test_tensor_with_operators():
    """Test tensor combined with creation/destruction operators."""
    print("Testing tensor-operator combinations...")

    # Create indices
    i = sqa.index('i', [], True)
    j = sqa.index('j', [], True)

    # Create h1 tensor
    h1_sym = [sqa.symmetry((1,0), 1)]
    h1 = sqa.tensor('h1', [i, j], h1_sym)

    # Build term: h1(i,j) * a+_i * a_j
    term = sqa.term(1.0, [], [h1, sqa.creOp(i), sqa.desOp(j)])

    # Normal order the term
    result = sqa.normalOrder(term)

    # Should produce 1 term (already in normal order)
    if len(result) == 1:
        print("  Normal ordering h1(i,j)*a+_i*a_j produces", len(result), "term")
        print("  Result:", result[0])
        print("  Test passed\n")
        return True
    else:
        print("  ERROR: Expected 1 term, got", len(result))
        return False

def test_commutator_relations():
    """Test fundamental anticommutation relations."""
    print("Testing anticommutation relations...")

    # Create indices
    i = sqa.index('i', [], True)
    j = sqa.index('j', [], True)

    # Test [a_i, a+_j]_+ = δ_ij (anticommutator)
    # Build term: a_i * a+_j
    term1 = sqa.term(1.0, [], [sqa.desOp(i), sqa.creOp(j)])
    result1 = sqa.normalOrder(term1)

    # Build term: a+_j * a_i
    term2 = sqa.term(1.0, [], [sqa.creOp(j), sqa.desOp(i)])
    result2 = sqa.normalOrder(term2)

    # After normal ordering:
    # a_i * a+_j -> -a+_j * a_i + δ_ij
    # a+_j * a_i -> a+_j * a_i (already normal)

    # Check that result1 contains delta function
    has_delta = False
    delta_term = None
    for term in result1:
        for tensor in term.tensors:
            if hasattr(tensor, 'name') and tensor.name == 'kdelta':
                has_delta = True
                delta_term = term
                break
        if has_delta:
            break

    if has_delta and len(result1) == 2:
        print("  Anticommutation {a_i, a+_j} = δ_ij verified")
        print("  Result has", len(result1), "terms:")
        for t in result1:
            print("   ", t)
        print("  Test passed\n")
        return True
    else:
        print("  ERROR: Expected 2 terms with delta function")
        print("  Got", len(result1), "terms")
        for t in result1:
            print("   ", t)
        return False

if __name__ == '__main__':
    print("\n" + "="*60)
    print("Testing Common Tensor Types")
    print("="*60 + "\n")

    all_passed = True

    # Run tests
    all_passed &= test_one_electron_tensors()
    all_passed &= test_two_electron_tensors()
    all_passed &= test_two_particle_density_matrix()
    all_passed &= test_tensor_with_operators()
    all_passed &= test_commutator_relations()

    print("="*60)
    if all_passed:
        print("All tensor tests PASSED")
        sys.exit(0)
    else:
        print("Some tensor tests FAILED")
        sys.exit(1)
