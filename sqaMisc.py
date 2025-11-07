#    file:  sqaMisc.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Defines some miscilaneous functions.
#
# (c) 2008-2009 Eric Neuscamman (eric.neuscamman@gmail.com)
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# In addition, any modification or use of this software should
# cite the following paper:
#
#   E. Neuscamman, T. Yanai, and G. K.-L. Chan.
#   J. Chem. Phys. 130, 124102 (2009)
#
# OPTIMIZATION (2025-11-06):
#   Added generateValidMatchings() function for
#   efficient k-matching generation via backtracking.
# AUTHOR: Konstantin Komarov (constlike@gmail.com).


from functools import cmp_to_key, lru_cache
from sqaOptions import options
import time

# Python 3 compatibility: cmp() function was removed in Python 3
def cmp(a, b):
    return (a > b) - (a < b)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


@lru_cache(maxsize=4096)
def _makeTuples_cached(n, inTuple):
  """Cached version that works with tuples."""
  inList = list(inTuple)
  outList = []
  if n == 0:
    return [tuple()]
  if n == len(inList):
    return [tuple(inList)]
  if n == 1:
    return [tuple([elem]) for elem in inList]
  tempList = []
  for i in range(len(inList)-n+1):
    tempList.append(inList[i])
  for i in range(len(tempList)):
    subList = _makeTuples_cached(n-1, tuple(inList[i+1:]))
    for sub in subList:
      outList.append(tuple([tempList[i]] + list(sub)))
  return outList

def makeTuples(n,inList):
  "Returns a list of all n-tuples that can be formed from the elements of inList.\n" + \
  "Warning:  behavior may be crazy if inList consists of mutable objects."
  # Convert to hashable types for caching
  try:
    inTuple = tuple(inList)
    result = _makeTuples_cached(n, inTuple)
    # Convert back to lists for compatibility
    return [list(t) for t in result]
  except TypeError:
    # Fallback for unhashable types (shouldn't happen often)
    outList = []
    if n == 0:
      return [[]]
    if n == len(inList):
      return [inList]
    if n == 1:
      for i in range(len(inList)):
        outList.append([inList[i]])
      return outList
    tempList = []
    for i in range(len(inList)-n+1):
      tempList.append(inList[i])
    for i in range(len(tempList)):
      subList = makeTuples(n-1,inList[i+1:])
      for j in range(len(subList)):
        outList.append([tempList[i]])
        for k in range(len(subList[j])):
          outList[-1].append(subList[j][k])
    return outList


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def generateValidMatchings(k, edges):
  """
  Generate k-matchings where no vertex is used more than once.

  Backtracking algorithm that generates only valid matchings.

  Args:
    k: Number of edges to select
    edges: List of (left_vertex, right_vertex) tuples

  Yields:
    Lists of k edges where each vertex appears at most once
  """
  def backtrack(start_idx, selected, used_left, used_right, remaining_k):
    if remaining_k == 0:
      yield list(selected)
      return

    for idx in range(start_idx, len(edges)):
      left, right = edges[idx]

      if left in used_left or right in used_right:
        continue

      selected.append(edges[idx])
      used_left.add(left)
      used_right.add(right)

      yield from backtrack(idx + 1, selected, used_left, used_right, remaining_k - 1)

      selected.pop()
      used_left.remove(left)
      used_right.remove(right)

  yield from backtrack(0, [], set(), set(), k)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def allDifferent(x):
  "Returns True if all the elements of x are different and False otherwise."
  # Optimized: O(n) instead of O(n²)
  return len(x) == len(set(x))


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def makePermutations(n):
  "Returns a list of all permutations of the order of the integers 0, 1, ..., n-1"
  if n == 1:
    return [[0]]
  outlist = []
  for i in range(n):
    temp = makePermutations(n-1)
    for j in range(len(temp)):
      for k in range(len(temp[j])):
        if temp[j][k] >= i:
          temp[j][k] += 1
      outlist.append(temp[j])
      outlist[-1].insert(0,i)
  return outlist


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def get_num_perms(ti,bi):
  """
  Returns the number of permutations between two lists of the integers 0 to N.
  """

  x = []
  for i in range(len(ti)):
    x.append([ti[i],bi[i]])
  x.sort(key=cmp_to_key(lambda a,b: cmp(a[1],b[1])))
  for i in range(len(x)):
    x[i] = x[i][0]
  i = 0
  n_perms = 0
  while i < len(x)-1:
    if x[i] != i:
      t = x[i]
      x[i] = x[t]
      x[t] = t
      n_perms += 1
    else:
      i += 1
  return n_perms


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------




def contractDeltaFuncs_complete(terms, zero_remaining=True):
    """
    Contract kdelta functions with option to zero remaining deltas.

    This is an extension of the built-in term.contractDeltaFuncs() method
    that optionally zeros out terms with remaining kdelta(i,j) where i!=j.

    Parameters
    ----------
    terms : list of term
        List of term objects to process
    zero_remaining : bool, default=True
        If True, zero out terms with remaining kdelta(i,j) where i!=j
        If False, keep terms with remaining kdelta (standard SQA behavior)

    Returns
    -------
    list of term
        Processed terms (original list is modified in-place)

    Notes
    -----
    The built-in contractDeltaFuncs() only removes diagonal kdelta(i,i)=1.
    Non-diagonal kdelta(i,j) where i!=j are left in the term.

    With zero_remaining=True, this function automatically zeros such terms,
    which is useful when evaluating overlap of orthogonal determinants.

    Examples
    --------
    >>> # Option 1: Zero remaining kdelta (default)
    >>> contractDeltaFuncs_complete(terms, zero_remaining=True)

    >>> # Option 2: Keep remaining kdelta (standard SQA behavior)
    >>> contractDeltaFuncs_complete(terms, zero_remaining=False)

    References
    ----------
    Added by Stan Pavlov (2025) for SI-SA-REKS verification.
    Extends SQA functionality without modifying core source code.
    """

    # Step 1: Call standard SQA contractDeltaFuncs for each term
    for term in terms:
        term.contractDeltaFuncs()

    # Step 2: If zero_remaining=True, zero out terms with remaining kdelta
    if zero_remaining:
        for term in terms:
            # Check if term has remaining kdelta functions
            has_kdelta = any(t.name == 'kdelta' for t in term.tensors)
            if has_kdelta:
                # Zero out this term (kdelta(i,j) with i!=j gives 0)
                term.numConstant = 0.0


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def filter_nonzero_terms(terms):
    """
    Filter out terms with zero coefficient.

    Parameters
    ----------
    terms : list of term
        List of term objects to filter

    Returns
    -------
    list of term
        Only non-zero terms

    Examples
    --------
    >>> final_terms = filter_nonzero_terms(norm_terms)

    References
    ----------
    Added by Stan Pavlov (2025) for SI-SA-REKS verification.
    """
    return [t for t in terms if abs(t.numConstant) > 1e-12]


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def filter_fully_contracted(terms):
    """
    Filter out terms that still have operators (keep only vacuum expectation values).

    Parameters
    ----------
    terms : list of term
        List of term objects to filter

    Returns
    -------
    list of term
        Only fully contracted terms (no operators)

    Examples
    --------
    >>> final_terms = filter_fully_contracted(normalized_terms)

    References
    ----------
    Added by Stan Pavlov (2025) for SI-SA-REKS verification.
    """
    from sqaTensor import creOp, desOp

    result = []
    for term in terms:
        has_ops = any(isinstance(t, (creOp, desOp)) for t in term.tensors)
        if not has_ops:
            result.append(term)
    return result


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

