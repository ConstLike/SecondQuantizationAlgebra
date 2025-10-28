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

## 6. Performance Notes

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
- We assume Mulliken notation for two-electron integrals

**Index ordering convention**:
- Creation operators ordered right-to-left
- Destruction operators ordered right-to-left
- Normal order: all creation operators left of destruction operators

---

*Second Quantization Algebra User Guide*

*Python 3 optimized version - October 2025*

*Maintained by Konstantin Komarov (constlike@gmail.com)*
