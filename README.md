# Second Quantization Algebra (SQA) Package

## Overview

The working equations of the full-internally contracted multireference configuration interaction (FIC-MRCI) were derived with the modified SQA program using spin-free unitary group generators. The code generator handles tensor contractions that cannot be coded manually. This feature enabled FIC-MRCI implementation as a pilot code for DMRG-MRCI development (M. Saitow et al. J Chem Phys 139, 044118 (2013)). The modified SQA package supports development of electronic structure theory methods.

## Current Version

**Maintainer**: Konstantin Komarov  
**Last Updated**: October 2025  
**Python Version**: 3.x (migrated and optimized from Python 2.x)

## Installation

### From Source

```bash
# Clone repository
cd /path/to/SecondQuantizationAlgebra

# Install in development mode
pip install -e .

# Or install directly
python setup.py install
```

### Quick Start

```python
import secondQuantizationAlgebra as sqa

# Create indices
i = sqa.index('i', [], True)
j = sqa.index('j', [], True)

# Create operators
cre_i = sqa.creOp(i)
des_j = sqa.desOp(j)

# Normal order
term = sqa.term(1.0, [], [des_j, cre_i])
result = sqa.normalOrder(term)
```

## Documentation

For detailed documentation with examples and API reference, see [SQA_User_Guide.md](SQA_User_Guide.md).

The User Guide covers:
- Core classes (indices, operators, tensors, terms)
- Operations (normal ordering, commutators, RDM decomposition)
- Code generation (Fortran 90, TeX output)
- Utility functions and workflows

## Testing

The package includes test suites to verify functionality:

### Run Basic Tests

```bash
python3 sqaTest.py
```

Tests 1-3 validate core operations:
- **Test 1**: Normal ordering of operator strings
- **Test 2**: Commutators with symmetry handling
- **Test 3**: Index type constraints

### Run Tensor Tests

```bash
python3 sqaTest_tensors.py
```

Validates common tensor types (h1, h2, a2, d1, d2) and symmetries.

### Run MRSF Benchmark

```bash
python3 sqaTest_mrsf.py
```

Mixed Reference Spin-Flip TDDFT: T11ab calculation (benchmark: 2.3 seconds, 22 terms).

## Version History

### Python 3 Migration and Performance Optimization (2025)

Python 3 migration and performance improvements:

  * **Python 3 compatibility**: Migrated from Python 2.x using 2to3 and manual fixes
  * **Performance**: 14.7x speedup over Python 2 implementation
  * **Code changes**:
    - Replaced `del()` loop operations with list comprehensions
    - Added LRU caching to frequently called functions
    - Reduced `allDifferent()` complexity from O(n^2) to O(n)
  * **Benchmark**: MRSF T11ab test executes in 2.3 s (was 33.6 s in Python 2)

### Original Modifications by M. Saitow

Distinctions from the original SQA package:

  * Use of spin-free unitary group generator is made fully functional
  * A simple code generator integrated inside SQA

## Recent Extensions (2025)

### Helper Functions for Delta Function Handling

Three new helper functions extend SQA functionality:

#### `contractDeltaFuncs_complete(terms, zero_remaining=True)`

Extension of built-in `term.contractDeltaFuncs()` that optionally zeros out terms with remaining kdelta(i,j) where i!=j.

**Parameters:**
- `terms`: list of term objects to process
- `zero_remaining`: if True (default), zero out terms with remaining non-diagonal kdelta

**Example:**
```python
import secondQuantizationAlgebra as sqa

# Calculate overlap <aa|bb> (should be zero)
# ... create terms with normalOrder() ...

# Option 1: Auto-zero remaining kdelta (default)
sqa.contractDeltaFuncs_complete(terms, zero_remaining=True)

# Option 2: Keep remaining kdelta (standard SQA behavior)
sqa.contractDeltaFuncs_complete(terms, zero_remaining=False)
```

**Note:** The built-in `contractDeltaFuncs()` only removes diagonal kdelta(i,i)=1. Non-diagonal kdelta(i,j) where i!=j are left in the term. This extension automatically zeros such terms, useful for evaluating overlap of orthogonal determinants.

#### `filter_nonzero_terms(terms)`

Filter out terms with zero coefficient (threshold: 1e-12).

```python
final_terms = sqa.filter_nonzero_terms(norm_terms)
```

#### `filter_fully_contracted(terms)`

Filter out terms that still have operators (keep only vacuum expectation values).

```python
vev_terms = sqa.filter_fully_contracted(normalized_terms)
```

---

## Citations

If you use or modify this version of SQA package in your own work, please cite the following publications:

### Original SQA Package

  * E. Neuscamman, T. Yanai, and G. K.-L. Chan, "Quadratic canonical transformation theory and higher order density matrices", *J. Chem. Phys.* **130**, 124102 (2009). DOI: 10.1063/1.3086932

### Spin-Free Extensions

  * M. Saitow, Y. Kurashige, and T. Yanai, "Multireference configuration interaction theory using cumulant reconstruction with internal contraction of density matrix renormalization group wave function", *J. Chem. Phys.* **139**, 044118 (2013). DOI: 10.1063/1.4816627

### Application Example

  * K. Komarov, H. Zhai, and G. K.-L. Chan, "Accurate Spin-Orbit Coupling by Relativistic Mixed-Reference Spin-Flip-TDDFT", *J. Chem. Theory Comput.* **19**(3), 953-964 (2023). DOI: 10.1021/acs.jctc.2c01036

---

## License

This software is provided under an open-source license. The original SQA package was developed by E. Neuscamman, T. Yanai, and G. K.-L. Chan. Modifications by M. Saitow and Python 3 migration by K. Komarov.

**Usage**: Free for academic and research purposes. Proper citation required.

**Warranty**: This program is distributed without warranty of any kind.

---

## Keywords

`second quantization`, `quantum chemistry`, `symbolic algebra`, `tensor algebra`, `coupled cluster`, `configuration interaction`, `density matrices`, `RDM`, `normal ordering`, `commutators`, `spin-free`, `multireference`, `MRCI`, `DMRG`, `Python 3`

