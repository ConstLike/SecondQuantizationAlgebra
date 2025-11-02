# Second Quantization Algebra (SQA) User Guide

## Overview

This guide covers the Second Quantization Algebra package. The package implements symbolic algebra for quantum chemistry using second quantization formalism to manipulate creation and destruction operators.

## Citations

If you use this package, please cite:

**Original SQA Package**:
- E. Neuscamman, T. Yanai, and G. K.-L. Chan, *J. Chem. Phys.* **130**, 124102 (2009). DOI: 10.1063/1.3086932

**Spin-Free Extensions**:
- M. Saitow, Y. Kurashige, and T. Yanai, *J. Chem. Phys.* **139**, 044118 (2013). DOI: 10.1063/1.4816627

**Application Example**:
- K. Komarov, H. Zhai, and G. K.-L. Chan, *J. Chem. Theory Comput.* **19**(3), 953-964 (2023). DOI: 10.1021/acs.jctc.2c01036

---

## 1. Core Classes

### 1.1 Index (`sqa.index`)

Define indices for tensor subscripts and operator labels.

**Syntax**: `sqa.index(name, indexType, isSummed, isExt)`

**Parameters**:
- `name` (str): Index label (e.g., 'i0', 'j1')
- `indexType` (list): Nested list specifying spin and orbital types
- `isSummed` (bool): True if index is summed (dummy index)
- `isExt` (bool): True if index represents external orbital

**Examples**:
```python
import secondQuantizationAlgebra as sqa

# Dummy index with alpha spin, all orbital types
i0 = sqa.index('i0', [['alpha'], ['core', 'active', 'virtual']], True, False)

# Fixed index with beta spin
j1 = sqa.index('j1', [['beta'], ['active']], False, False)

# Simple index without type constraints
k = sqa.index('k', [], True, False)
```

---

### 1.2 Operators

Fermionic creation and destruction operators.

#### Creation Operator (`sqa.creOp`)

**Syntax**: `sqa.creOp(index)`

**Example**:
```python
i = sqa.index('i', [], True)
cre_i = sqa.creOp(i)  # Creates a^+_i
```

#### Destruction Operator (`sqa.desOp`)

**Syntax**: `sqa.desOp(index)`

**Example**:
```python
j = sqa.index('j', [], True)
des_j = sqa.desOp(j)  # Creates a_j
```

---

### 1.3 Tensors (`sqa.tensor`)

We define general tensors with symmetry properties.

**Syntax**: `sqa.tensor(name, indices, symmetries)`

**Parameters**:
- `name` (str): Tensor identifier
- `indices` (list): List of index objects
- `symmetries` (list): List of symmetry objects

**Common tensors**:
- `h1`: One-electron integrals
- `h2`: Two-electron integrals
- `a2`: Two-electron amplitudes
- `d1`: One-particle density matrix
- `d2`: Two-particle density matrix

**Example**:
```python
# One-electron Hamiltonian with exchange symmetry
i0 = sqa.index('i0', [], True)
i1 = sqa.index('i1', [], True)
h1_sym = [sqa.symmetry((1,0), 1)]  # h(i,j) = h(j,i)
h1 = sqa.tensor('h1', [i0, i1], h1_sym)
```

#### Two-Electron Integrals and Mulliken Notation

This package uses Mulliken notation for two-electron integrals. This is important convention that affects how integrals combine with creation and destruction operators.

**Mulliken notation definition**:

The two-electron integral in Mulliken notation (ij|kl) represents:

```
(ij|kl) = integral integral phi_i(r1) phi_j(r1) (1/r12) phi_k(r2) phi_l(r2) dr1 dr2
```

where r12 is distance between electrons 1 and 2.

**Symmetry properties**:

Mulliken integrals have these symmetries:
- (ij|kl) = (ji|kl) for non-antisymmetrized integrals
- (ij|kl) = (kl|ij) (particle exchange symmetry)
- (ij|kl) = (ij|lk) for non-antisymmetrized integrals

For antisymmetrized integrals used in this package:
- (ij|kl) = -(ji|kl) (antisymmetric in first pair)
- (ij|kl) = -(ij|lk) (antisymmetric in second pair)
- (ij|kl) = (kl|ij) (particle exchange)

**Critical convention - operator ordering**:

When we construct Hamiltonian terms with two-electron integrals, destruction operators appear in **reverse order** relative to last two integral indices:

```python
# For integral h2(i,j,k,l) representing (ij|kl)
# Operators are: a+_i a+_j a_l a_k (note: a_l a_k, not a_k a_l)
```

This is **mandatory** convention in this package.

**Example - antisymmetrized two-electron integrals**:

```python
# Create indices
p0 = sqa.index('p0', [], True)
p1 = sqa.index('p1', [], True)
p2 = sqa.index('p2', [], True)
p3 = sqa.index('p3', [], True)

# Define symmetries for antisymmetrized Mulliken integrals
# Antisymmetry in first index pair (indices 0,1)
sym_antisym_12 = sqa.symmetry((1,0,2,3), -1)  # (ij|kl) = -(ji|kl)
# Exchange symmetry between index pairs (0,1) <-> (2,3)
sym_exch_pairs = sqa.symmetry((2,3,0,1), 1)   # (ij|kl) = (kl|ij)

# Create h2 tensor
h2 = sqa.tensor('h2', [p0, p1, p2, p3], [sym_antisym_12, sym_exch_pairs])

# Build Hamiltonian term
# IMPORTANT: destruction operators in reverse order!
# h2(p0,p1,p2,p3) corresponds to integral (p0 p1|p2 p3)
# Operators: a+_p0 a+_p1 a_p3 a_p2 (note a_p3 a_p2 order)
H_term = sqa.term(1.0, [], [h2,
                             sqa.creOp(p0),
                             sqa.creOp(p1),
                             sqa.desOp(p3),  # p3 comes before p2
                             sqa.desOp(p2)])
```

**Why this convention**:

This ordering follows from second quantization theory where two-electron operator is:

```
H = (1/2) sum_ijkl (ij|kl) a+_i a+_j a_l a_k
```

The destruction operator order (a_l a_k instead of a_k a_l) comes from normal ordering conventions and ensures correct signs when evaluating matrix elements.

**Common mistake**:

```python
# WRONG - destruction operators in wrong order
H_wrong = sqa.term(1.0, [], [h2,
                              sqa.creOp(p0),
                              sqa.creOp(p1),
                              sqa.desOp(p2),  # Wrong!
                              sqa.desOp(p3)])

# CORRECT - destruction operators reversed
H_correct = sqa.term(1.0, [], [h2,
                                sqa.creOp(p0),
                                sqa.creOp(p1),
                                sqa.desOp(p3),  # Correct
                                sqa.desOp(p2)])
```

**Verification**:

You can verify Mulliken convention by checking that your integrals satisfy:
- h2(i,j,k,l) = -h2(j,i,k,l) for antisymmetrized case
- h2(i,j,k,l) = h2(k,l,i,j) always

These symmetries are implemented through symmetry objects shown in example above.

---

### 1.4 Symmetries (`sqa.symmetry`)

We specify tensor symmetry under index permutations.

**Syntax**: `sqa.symmetry(permutation, sign)`

**Parameters**:
- `permutation` (tuple): New index order (0-indexed)
- `sign` (int): +1 for symmetric, -1 for antisymmetric

**Examples**:
```python
# Symmetric exchange: T(i,j) = T(j,i)
sym1 = sqa.symmetry((1,0), 1)

# Antisymmetric exchange: T(i,j,k,l) = -T(j,i,k,l)
sym2 = sqa.symmetry((1,0,2,3), -1)

# Particle exchange: T(i,j,k,l) = T(k,l,i,j)
sym3 = sqa.symmetry((2,3,0,1), 1)
```

---

### 1.5 Terms (`sqa.term`)

We construct terms as products of tensors and operators.

**Syntax**: `sqa.term(coefficient, constants, tensors)`

**Parameters**:
- `coefficient` (float): Numerical prefactor
- `constants` (list): List of constant tensors
- `tensors` (list): List of tensor and operator objects

**Example**:
```python
# Construct term: 0.5 * h1(i,j) * a^+_i * a_j
i = sqa.index('i', [], True)
j = sqa.index('j', [], True)
h1 = sqa.tensor('h1', [i,j], [])
term1 = sqa.term(0.5, [], [h1, sqa.creOp(i), sqa.desOp(j)])
```

---

## 2. Core Operations

### 2.1 Normal Ordering

We transform operator strings to normal-ordered form.

**Function**: `sqa.normalOrder(term)`

**Returns**: List of terms in normal-ordered form

**Example**:
```python
# Normal order: a_i a_j a^+_k a^+_l
i = sqa.index('i', [], False)
j = sqa.index('j', [], False)
k = sqa.index('k', [], False)
l = sqa.index('l', [], False)

term_input = sqa.term(1.0, [], [
    sqa.desOp(i),
    sqa.desOp(j),
    sqa.creOp(k),
    sqa.creOp(l)
])

result = sqa.normalOrder(term_input)
# Returns 7 terms including contractions
```

**Maximum Contraction**: `sqa.normalOrder_maxcontraction(term)`

We compute only maximally contracted terms.

---

### 2.2 Commutators

We evaluate commutators between operators.

**Function**: `sqa.commutator(term1, terms2)`

**Parameters**:
- `term1` (term): Single term object
- `terms2` (list): List of term objects

**Returns**: List of terms representing [term1, sum(terms2)]

**Example**:
```python
# Compute [H, T2] where H = h1 a^+ a, T2 = t2 a^+ a^+ a a
i = sqa.index('i', [], True)
j = sqa.index('j', [], True)

h1 = sqa.tensor('h1', [i,j], [])
h_term = sqa.term(1.0, [], [h1, sqa.creOp(i), sqa.desOp(j)])

# Define T2 terms (two terms due to antisymmetry)
k = sqa.index('k', [], True)
l = sqa.index('l', [], True)
m = sqa.index('m', [], True)
n = sqa.index('n', [], True)
t2 = sqa.tensor('t2', [k,l,m,n], [])

t2_terms = []
t2_terms.append(sqa.term(1.0, [], [t2, sqa.creOp(k), sqa.creOp(l),
                                     sqa.desOp(n), sqa.desOp(m)]))
t2_terms.append(sqa.term(-1.0, [], [t2, sqa.creOp(m), sqa.creOp(n),
                                      sqa.desOp(l), sqa.desOp(k)]))

comm_result = sqa.commutator(h_term, t2_terms)
```

---

### 2.3 Term Manipulation

#### Combine Like Terms

**Function**: `sqa.combineTerms(terms)`

We merge terms with identical tensor structure.

**Example**:
```python
terms = sqa.normalOrder(some_term)
sqa.combineTerms(terms)  # Modifies list in-place
```

#### Multiply Terms

**Function**: `sqa.multiplyTerms(term1, term2, numNewIndices)`

We compute the product of two terms.

**Example**:
```python
result = sqa.multiplyTerms(term1, term2, 10)
```

#### Remove Small Terms

**Function**: `sqa.termChop(terms, threshold)`

We eliminate terms below numerical threshold.

**Example**:
```python
sqa.termChop(terms, 1e-10)  # Remove |coeff| < 1e-10
```

---

## 3. Advanced Features

### 3.1 Spin-Free Excitation Operators

We use unitary group generators for spin-free formulations.

**Class**: `sqa.sfExOp`

**Syntax**: `sqa.sfExOp(indices)`

**Example**:
```python
# Create E^{ij}_{kl} operator
indices = [sqa.index(s, [], True) for s in ['i','j','k','l']]
E_op = sqa.sfExOp(indices)
```

---

### 3.2 Operator to RDM Conversion

**Function**: `sqa.convert_ops_to_rdms_so(terms, rdm_name, ord=0)`

Converts normal ordered creation and destruction operators to corresponding reduced density matrices. Terms must be normal ordered before calling this function.

**Parameters**:
- `terms`: List of term objects with normal ordered operators
- `rdm_name`: String name for the RDM tensors (e.g., 'd1', 'd2')
- `ord`: If positive, only converts RDMs of that order (0 converts all)

**Example**:
```python
# Create terms with creation/destruction operators
i = sqa.index('i', [], True)
j = sqa.index('j', [], True)
term = sqa.term(1.0, [], [sqa.creOp(i), sqa.desOp(j)])

# Normal order first
terms = sqa.normalOrder(term)

# Convert operators to RDM
sqa.convert_ops_to_rdms_so(terms, 'd1')
# Now terms contain d1(i,j) instead of a+_i a_j
```

---

### 3.3 RDM Decomposition

Higher-order reduced density matrices (RDMs) can be decomposed into products of lower-order RDMs. The n-particle RDM is defined as:

```
n-RDM(i1...in; j1...jn) = <Psi| a+_i1 ... a+_in a_jn ... a_j1 |Psi>
```

Decomposition uses cumulant expansion theory. The 3-RDM decomposes into products of 2-RDM and 1-RDM. The 4-RDM decomposes into 3-RDM, 2-RDM, and 1-RDM products. This reduces computational cost by expressing higher-order objects through lower-order ones.

#### Spin-Free Decomposition

**Functions**:
- `sqa.decomp_3rdm_to_2rdm_sf(terms)`: 3-RDM -> 2-RDM
- `sqa.decomp_4rdm_to_2rdm_sf(terms)`: 4-RDM -> 2-RDM
- `sqa.decomp_4rdm_to_3rdm_sf(terms)`: 4-RDM -> 3-RDM

#### Spin-Orbital Decomposition

**Functions**:
- `sqa.decomp_3rdm_to_2rdm_so(terms)`: 3-RDM -> 2-RDM
- `sqa.decomp_4rdm_to_2rdm_so(terms)`: 4-RDM -> 2-RDM

**Example**:
```python
# Decompose 3-RDM terms into 2-RDM contractions
terms_3rdm = [...]  # List of terms with 3-body RDMs
terms_2rdm = sqa.decomp_3rdm_to_2rdm_sf(terms_3rdm)
```

---

### 3.4 Code Generation

We convert symbolic expressions to executable code.

#### Fortran 90 Output

**Functions**:
- `sqa.convertF90_0(terms, filename)`: Zero-body terms
- `sqa.convertF90_1(terms, filename)`: One-body terms
- `sqa.convertF90_ERI(terms, filename)`: Two-body terms

#### TeX Output

**Function**: `sqa.convertTeX(terms, filename)`

We generate LaTeX-formatted equations.

**Example**:
```python
terms = sqa.commutator(h_term, t2_terms)
sqa.combineTerms(terms)
sqa.convertTeX(terms, 'output.tex')
```

---

## 4. Utility Functions

### 4.1 Index and Tensor Utilities

**Transpose Equivalence**: `sqa.combine_transpose(terms)`

We identify and merge transpose-equivalent terms.

**RDM Type Assignment**: `sqa.assign_rdm_types(term, rdm_name)`

We assign orbital types to RDM indices.

---

### 4.2 List Operations

**Make Tuples**: `sqa.makeTuples(n, inList)`

We generate all n-element combinations from list.

**Example**:
```python
result = sqa.makeTuples(2, [1,2,3,4])
# Returns: [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]]
```

**Make Permutations**: `sqa.makePermutations(n)`

We generate all permutations of integers 0 to n-1.

**All Different**: `sqa.allDifferent(x)`

We check if all list elements are unique.

---

### 4.3 Operator Filtering

**Sort Operators**: `sqa.sortOps(unsortedOps, returnPermutation=False)`

Sorts creation and destruction operators into normal order without performing contractions.

**Example**:
```python
ops = [sqa.desOp(i), sqa.creOp(j), sqa.desOp(k)]
(sign, sortedOps) = sqa.sortOps(ops)
# Returns sign from permutation and sorted operator list
```

**Remove Core Operator Pairs**: `sqa.removeCoreOpPairs(termList)`

Removes pairs of core creation and destruction operators for the same index. Terms must be normal ordered before calling this function.

**Example**:
```python
terms = sqa.normalOrder(input_term)
sqa.removeCoreOpPairs(terms)
```

**Remove Core Operators (Spin-Free)**: `sqa.removeCoreOps_sf(termList)`

Removes core indices from spin-free excitation operators. Use this before converting operators to density matrices. Terms with zero expectation values are removed.

**Remove Virtual Operators (Spin-Free)**: `sqa.removeVirtOps_sf(termList)`

Removes terms containing spin-free operators with virtual indices.

---

### 4.4 Function Dependencies

Several functions call other functions automatically. Understanding these dependencies prevents redundant calls:

**combineTerms automatically calls**:
- `makeCanonical()` on each term (converts to canonical index naming)
- `sort()` on the term list (orders terms for comparison)
- `termChop()` at the end (removes zero coefficients)

You should NOT call these functions separately after using `combineTerms`.

**normalOrder does NOT automatically call**:
- `combineTerms()` - you must call this explicitly to merge like terms

**Typical workflow**:
```python
# Step 1: Normal order
terms = sqa.normalOrder(input_term)

# Step 2: Combine like terms (calls makeCanonical and termChop automatically)
sqa.combineTerms(terms)

# Step 3: Apply any filtering if needed
sqa.removeCoreOpPairs(terms)
```

**Code generation workflow**:
```python
# Derive equations
terms = sqa.commutator(h_term, t_terms)
sqa.combineTerms(terms)  # Combines and canonicalizes automatically

# Remove unnecessary operators
sqa.removeCoreOpPairs(terms)

# Generate code
sqa.convertF90_1(terms, 'output.f90')
```

---

## 5. Common Workflows

### 5.1 Derive CCSD Equations

```python
import secondQuantizationAlgebra as sqa

# Define indices
indices = [sqa.index('i%d'%n, [], True) for n in range(10)]

# Define Hamiltonian
h1 = sqa.tensor('f', indices[0:2], [])
h2 = sqa.tensor('v', indices[2:6], [])

# Define cluster operators T1 and T2
t1 = sqa.tensor('t1', indices[0:2], [])
t2 = sqa.tensor('t2', indices[0:4], [])

# Build operator strings
H_term = [...]  # Hamiltonian terms
T1_terms = [...]  # T1 cluster operator
T2_terms = [...]  # T2 cluster operator

# Compute [H, T1] and [H, T2]
comm1 = sqa.commutator(H_term, T1_terms)
comm2 = sqa.commutator(H_term, T2_terms)

# Combine and simplify
all_terms = comm1 + comm2
sqa.combineTerms(all_terms)
```

### 5.2 Normal Order Product

```python
# Create operator product
i, j, k, l = [sqa.index(s, [], False) for s in ['i','j','k','l']]
ops = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(k), sqa.desOp(l)]
term = sqa.term(1.0, [], ops)

# Normal order with all contractions
result = sqa.normalOrder(term)
sqa.combineTerms(result)

# Print result
for t in result:
    print(t)
```

### 5.3 Multireference Calculations

```python
# Use spin-free generators
E_ops = []
for i in range(n_states):
    indices = [...]  # Define state indices
    E_ops.append(sqa.sfExOp(indices))

# Contract with Hamiltonian
H_terms = [...]
result = []
for E in E_ops:
    term = sqa.term(1.0, [], [E] + H_terms)
    contracted = sqa.contractCoreOps_sf(term)
    result.extend(contracted)
```

---

## 6. Parallel Computation with sqaParallel

For computationally intensive matrix element calculations, the `sqaParallel` module provides automatic parallelization using Python's multiprocessing. This bypasses the Global Interpreter Lock (GIL) to achieve parallel speedup on multi-core systems.

### 6.1 Overview

The sqaParallel module is designed for computing matrix elements between configurations (e.g., Slater determinants) in parallel. It provides:

- **Automatic CPU detection**: Uses all available cores by default
- **Hermitian symmetry**: Exploits <I|H|J> = <J|H|I> symmetry (1.6x reduction)
- **Selection rules**: Optional Slater-Condon rules to skip zero elements
- **Progress tracking**: Real-time reporting with process IDs
- **Result verification**: Built-in symmetry checking

**Typical speedups**: 4-8x on 4-10 core systems, depending on matrix size and optimizations enabled.

### 6.2 Basic Usage

The core function is `compute_matrix()`, which takes a user-defined computation function and parallelizes it across all matrix elements.

**Example**:

```python
import secondQuantizationAlgebra as sqa
import sqaParallel

# Define your computation function (MUST be at module level for pickling)
def compute_element(bra_name, ket_name):
    """
    Compute <bra|H|ket> matrix element.

    This function is called in parallel for each matrix element.
    Returns: list of SQA term objects
    """
    # Your SQA computation here
    # 1. Build bra and ket as SQA terms
    # 2. Build Hamiltonian operator
    # 3. Compute matrix element using normalOrder()
    # 4. Contract and combine terms

    result_terms = [...]  # list of sqa.term objects
    return result_terms

# Define configurations
config_names = ['config1', 'config2', 'config3', 'config4']

# Compute matrix with automatic parallelization
matrix, stats = sqaParallel.compute_matrix(
    compute_func=compute_element,
    config_names=config_names,
    use_symmetry=True,        # Use Hermitian symmetry
    selection_rule_func=None, # Optional selection rules
    n_workers=None,           # Auto-detect CPU count
    verbose=True              # Print progress
)

# Access results
result = matrix[('config1', 'config2')]  # List of terms
```

**Important**: The `compute_func` must be defined at **module level** (not inside another function or class) for Python's multiprocessing to serialize it correctly.

### 6.3 Configuration Options

**Parameters for `compute_matrix()`**:

- **`compute_func`**: User-defined function with signature `result = func(bra_name, ket_name)`
  - Must return list of SQA term objects
  - Must be defined at module level for pickling

- **`config_names`**: List of configuration names (e.g., `['aabb', 'aacc', ...]`)

- **`use_symmetry`**: bool, default=True
  - Exploit Hermitian symmetry: compute upper triangle + diagonal only
  - Copy lower triangle from upper triangle
  - Reduces computation by ~1.6x

- **`selection_rule_func`**: callable or None, default=None
  - Optional function to skip zero elements by selection rules
  - Signature: `is_zero = func(bra_name, ket_name)`
  - Returns True if element is known to be zero

- **`n_workers`**: int or None, default=None
  - Number of parallel workers
  - If None, auto-detect using `cpu_count()`
  - If 0, run sequentially (for debugging)

- **`verbose`**: bool, default=True
  - Print progress messages with process IDs

**Returns**:
- **`matrix`**: dict with (bra_name, ket_name) keys and term lists as values
- **`statistics`**: dict with timing info and speedup estimates

### 6.4 Selection Rules

The module provides helper functions for Slater-Condon selection rules, which state that matrix elements are zero if configurations differ by too many orbitals.

#### Slater-Condon Difference

**Function**: `sqaParallel.slater_condon_diff(config_a, config_b)`

Computes the number of spin-orbitals that differ between two configurations:

```
diff = |A \ B| = |B \ A|
```

where A and B are sets of occupied spin-orbitals.

**Example**:
```python
# Configuration aabb: ['a_alpha', 'a_beta', 'b_alpha', 'b_beta']
# Configuration aacc: ['a_alpha', 'a_beta', 'c_alpha', 'c_beta']

diff = sqaParallel.slater_condon_diff(
    ['a_alpha', 'a_beta', 'b_alpha', 'b_beta'],
    ['a_alpha', 'a_beta', 'c_alpha', 'c_beta']
)
# Returns: 2 (b_alpha, b_beta in first; c_alpha, c_beta in second)
```

#### Two-Electron Selection Rules

**Function**: `sqaParallel.make_selection_rule_two_electron(configurations, threshold=3)`

Creates a selection rule function for two-electron operators (e.g., H2e).

For two-electron operators:
```
<I|H2e|J> = 0  if  diff(I,J) >= 3
```

**Example**:
```python
# Define configurations
configurations = {
    'aabb': ['a_alpha', 'a_beta', 'b_alpha', 'b_beta'],
    'aacc': ['a_alpha', 'a_beta', 'c_alpha', 'c_beta'],
    'bbdd': ['b_alpha', 'b_beta', 'd_alpha', 'd_beta'],
    'ccdd': ['c_alpha', 'c_beta', 'd_alpha', 'd_beta']
}

# Create selection rule function
selection_rule = sqaParallel.make_selection_rule_two_electron(
    configurations,
    threshold=3
)

# Use in compute_matrix
matrix, stats = sqaParallel.compute_matrix(
    compute_func=compute_element,
    config_names=list(configurations.keys()),
    use_symmetry=True,
    selection_rule_func=selection_rule,  # Skip elements with diff >= 3
    n_workers=None
)
```

**Custom selection rules**: You can define your own selection rule function:

```python
def my_selection_rule(bra_name, ket_name):
    """Return True if element should be skipped (is zero)."""
    # Custom logic here
    if some_condition(bra_name, ket_name):
        return True  # Skip this element
    return False     # Compute this element
```

### 6.5 Performance Characteristics

**Granularity**: The sqaParallel module parallelizes at the matrix element level. Each worker computes complete matrix elements (which may take 10-30 seconds each for complex operators).

Do **NOT** try to parallelize individual `normalOrder()` calls, as:
- Each normalOrder() call takes 3-20 ms
- Multiprocessing overhead is 10-50 ms per task
- This would cause 10-100x slowdown instead of speedup

**Expected speedups** (4-electron REKS(4,4) example):

| Optimization | Speedup | Time (4x4 matrix) |
|--------------|---------|-------------------|
| None (sequential) | 1.0x | ~26 min |
| Symmetry only | 1.6x | ~16 min |
| Symmetry + Selection | 2.0x | ~13 min |
| Symmetry + Selection + Parallel (10 cores) | 7.9x | ~3.3 min |

**Efficiency**: Near-linear scaling up to ~10 cores for typical matrix sizes (4x4 to 16x16).

### 6.6 Verification

The module includes a verification function to check Hermitian symmetry:

```python
# Verify all elements satisfy <I|H|J> = <J|H|I>
is_symmetric = sqaParallel.verify_symmetry(
    matrix,
    config_names,
    verbose=True
)
```

### 6.7 Complete Example

For a complete working example, see `sqaTest_multiprocessing.py`:

```bash
python3 sqaTest_multiprocessing.py
```

This test computes diagonal matrix elements for a 4-electron system and demonstrates:
- Sequential vs parallel execution
- Result verification
- Speedup measurement

Expected output: ~4x speedup on 4-core systems (30s parallel vs 120s sequential).

---

## 7. Performance Notes

Key functions optimized for Python 3 compatibility:

- **Normal ordering**: 14.7x faster than original Python 2 version
- **Caching**: LRU cache applied to tuple generation
- **List operations**: Optimized comprehensions replace loop deletions

**Benchmark**: MRSF T11ab test executes in 2.3 s (was 33.6 s in Python 2)

---

## 7. Limitations and Notes

- We support only fermionic operators (no bosons)
- We require antisymmetrized integrals for proper results
- We use Einstein summation convention for repeated indices
- We assume Mulliken notation for two-electron integrals (see Section 1.3 for details)

**Index ordering convention**:
- Creation operators ordered right-to-left
- Destruction operators ordered right-to-left
- Normal order: all creation operators left of destruction operators
- Two-electron integrals: destruction operators in **reverse** order (see Section 1.3)

---

*Second Quantization Algebra User Guide*

*Python 3 optimized version - October 2025*

*Maintained by Konstantin Komarov (constlike@gmail.com)*
