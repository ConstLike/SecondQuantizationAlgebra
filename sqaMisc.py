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


