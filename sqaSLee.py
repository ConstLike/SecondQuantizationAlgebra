#    file:  sqaSLee.py
#  author:  Seunhoon Lee 
#    date:  Feb 30, 2021
# summary:  Various methods to get (T) correction 
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


from sqaTensor import tensor, kroneckerDelta, creOp, desOp, sfExOp
from sqaTerm import term, sortOps
from sqaIndex import index
from sqaMisc import makeTuples, allDifferent, makePermutations, get_num_perms
from sqaOptions import options


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

def HFFermi(inTerm):
  ""

#  if options.verbose:
#    print "converting to normal order:  %s" %(str(inTerm))

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError, "inTerm must be of class term"

  # determine what types of operators the term contains
  has_creDesOps = False
  has_sfExOps = False
  for t in inTerm.tensors:
    if isinstance(t, creOp) or isinstance(t, desOp):
      has_creDesOps = True
    elif isinstance(t, sfExOp):
      has_sfExOps = True

  # If not spin free excitation operators,
  # raise an error
  if not has_sfExOps:
    raise RuntimeError, "Present version of HFFermi can treat only sfExOp"

  if has_creDesOps and has_sfExOps:
    raise RuntimeError, "Present version of HFFermi can treat only sfExOp"

  # if the term is already normal ordered, return it unchanged
  elif not inTerm.isNormalOrdered():
    raise RuntimeError, "HFFermi requires normal ordered operators"

#  # Normal ordering for creOp/desOp
#  elif has_creDesOps:
#
#    # Separate the cre/des operators from other tensors
#    ops = []
#    nonOps = []
#    for t in inTerm.tensors:
#      if isinstance(t, creOp) or isinstance(t, desOp):
#        ops.append(t.copy())
#      else:
#        nonOps.append(t.copy())
#
#    # Generate all contraction pairs
#    contractionPairs = []
#    for i in range(len(ops)):
#      iTerm = ops[i]
#      for j in range(i+1,len(ops)):
#        jTerm = ops[j]
#        if isinstance(iTerm, desOp) and isinstance(jTerm, creOp):
#          contractionPairs.append((i,j));
#    #print "contractionPairs\n", contractionPairs
#
#    # Determine maximum contraction order
#    creCount = 0
#    maxConOrder = 0
#    for i in range(len(ops)-1,-1,-1):
#      iTerm = ops[i]
#      if isinstance(iTerm, creOp):
#        creCount +=1
#      elif isinstance(iTerm, desOp) and creCount > 0:
#        maxConOrder += 1
#        creCount -= 1
#    del(creCount,iTerm)
#
#    # Generate all contractions
#    contractions = []
#    for i in range(maxConOrder+1):
#      subCons = makeTuples(i,contractionPairs)
#      j = 0
#      while j < len(subCons):
#        creOpTags = []
#        desOpTags = []
#        for k in range(i):
#          creOpTags.append(subCons[j][k][1])
#          desOpTags.append(subCons[j][k][0])
#        if allDifferent(creOpTags) and allDifferent(desOpTags):
#          j += 1
#        else:
#          del(subCons[j])
#      for j in range(len(subCons)):
#        contractions.append(subCons[j])
#    del(subCons,creOpTags,desOpTags,contractionPairs)
#    #print "contractions:\n", contractions
#
#    # For each contraction, generate the resulting term
#    outTerms = []
#    for contraction in contractions:
#      conSign = 1
#      deltaFuncs = []
#      conIndeces = []
#      subOpString = []
#      subOpString.extend(ops)
#      for conPair in contraction:
#        index1 = ops[conPair[0]].indices[0]
#        index2 = ops[conPair[1]].indices[0]
#        deltaFuncs.append(kroneckerDelta([index1,index2]))
#        subOpString[conPair[0]] = 'contracted'
#        subOpString[conPair[1]] = 'contracted'
#        for q in subOpString[conPair[0]+1:conPair[1]]:
#          if not (q is 'contracted'):
#            conSign *= -1
#      i = 0
#      while i < len(subOpString):
#        if subOpString[i] is 'contracted':
#          del(subOpString[i])
#        else:
#          i += 1
#      (sortSign,sortedOps) = sortOps(subOpString)
#      totalSign = conSign * sortSign
#      outTensors = []
#      outTensors.extend(nonOps)
#      outTensors.extend(deltaFuncs)
#      outTensors.extend(sortedOps)
#      outTerms.append( term(totalSign * inTerm.numConstant, inTerm.constants, outTensors) )

  # Normal ordering for sfExOps
  elif has_sfExOps:

    # Make separate lists of the spin free excitation operators and other tensors
    sfExOp_list = []
    other_list  = []
    for t in inTerm.tensors:
      if isinstance(t, sfExOp):
        sfExOp_list.append(t.copy())
      else:
        other_list.append(t.copy())

    # Initialize n, the number of remaining spin free excitation operators
    n = len(sfExOp_list)
    if n != 1:
        raise RuntimeError, "terms should have single sfExOp"

    # 
    const = inTerm.numConstant
    order = sfExOp_list[0].order
    for i in range(order): 
        idx = sfExOp_list[0].indices[i+order] 
        #print idx.indType
        for typ in idx.indType:
            #print typ[0]
            if typ == options.virtual_type:
                const = 0.0 

    deltaFuncs = []
    if const == inTerm.numConstant:
        for i in range(order):
            index1 = sfExOp_list[0].indices[i]
            index2 = sfExOp_list[0].indices[i+order]
            deltaFuncs.append(kroneckerDelta([index1,index2])) 
            const *= 2
        outTerm = term(const, inTerm.constants, other_list + deltaFuncs)
    else:
        outTerm = term(const, inTerm.constants, other_list + sfExOp_list)

  else:
    raise RuntimeError, "HFFermi failed to choose what to do."

  return outTerm

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

def idxname(inTerm):
  """
  1. change the index form (core: i_k, active: r_k, virtual: a_k, whole: m_k)
  """ 
  

#  if options.verbose:
#    print "converting to normal order:  %s" %(str(inTerm))

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError, "inTerm must be of class term"

#  # determine what types of operators the term contains
#  has_creDesOps = False
#  has_sfExOps = False
#  for t in inTerm.tensors:
#    if isinstance(t, creOp) or isinstance(t, desOp):
#      has_creDesOps = True
#    elif isinstance(t, sfExOp):
#      has_sfExOps = True
#
#  # If not spin free excitation operators,
#  # raise an error
#  if not has_sfExOps:
#    raise RuntimeError, "Present version of HFFermi can treat only sfExOp"
#
#  if has_creDesOps and has_sfExOps:
#    raise RuntimeError, "Present version of HFFermi can treat only sfExOp"

  # if the term is already normal ordered, return it unchanged
  elif not inTerm.isNormalOrdered():
    raise RuntimeError, "HFFermi requires normal ordered operators"

#  # Normal ordering for creOp/desOp
#  elif has_creDesOps:
#
#    # Separate the cre/des operators from other tensors
#    ops = []
#    nonOps = []
#    for t in inTerm.tensors:
#      if isinstance(t, creOp) or isinstance(t, desOp):
#        ops.append(t.copy())
#      else:
#        nonOps.append(t.copy())
#
#    # Generate all contraction pairs
#    contractionPairs = []
#    for i in range(len(ops)):
#      iTerm = ops[i]
#      for j in range(i+1,len(ops)):
#        jTerm = ops[j]
#        if isinstance(iTerm, desOp) and isinstance(jTerm, creOp):
#          contractionPairs.append((i,j));
#    #print "contractionPairs\n", contractionPairs
#
#    # Determine maximum contraction order
#    creCount = 0
#    maxConOrder = 0
#    for i in range(len(ops)-1,-1,-1):
#      iTerm = ops[i]
#      if isinstance(iTerm, creOp):
#        creCount +=1
#      elif isinstance(iTerm, desOp) and creCount > 0:
#        maxConOrder += 1
#        creCount -= 1
#    del(creCount,iTerm)
#
#    # Generate all contractions
#    contractions = []
#    for i in range(maxConOrder+1):
#      subCons = makeTuples(i,contractionPairs)
#      j = 0
#      while j < len(subCons):
#        creOpTags = []
#        desOpTags = []
#        for k in range(i):
#          creOpTags.append(subCons[j][k][1])
#          desOpTags.append(subCons[j][k][0])
#        if allDifferent(creOpTags) and allDifferent(desOpTags):
#          j += 1
#        else:
#          del(subCons[j])
#      for j in range(len(subCons)):
#        contractions.append(subCons[j])
#    del(subCons,creOpTags,desOpTags,contractionPairs)
#    #print "contractions:\n", contractions
#
#    # For each contraction, generate the resulting term
#    outTerms = []
#    for contraction in contractions:
#      conSign = 1
#      deltaFuncs = []
#      conIndeces = []
#      subOpString = []
#      subOpString.extend(ops)
#      for conPair in contraction:
#        index1 = ops[conPair[0]].indices[0]
#        index2 = ops[conPair[1]].indices[0]
#        deltaFuncs.append(kroneckerDelta([index1,index2]))
#        subOpString[conPair[0]] = 'contracted'
#        subOpString[conPair[1]] = 'contracted'
#        for q in subOpString[conPair[0]+1:conPair[1]]:
#          if not (q is 'contracted'):
#            conSign *= -1
#      i = 0
#      while i < len(subOpString):
#        if subOpString[i] is 'contracted':
#          del(subOpString[i])
#        else:
#          i += 1
#      (sortSign,sortedOps) = sortOps(subOpString)
#      totalSign = conSign * sortSign
#      outTensors = []
#      outTensors.extend(nonOps)
#      outTensors.extend(deltaFuncs)
#      outTensors.extend(sortedOps)
#      outTerms.append( term(totalSign * inTerm.numConstant, inTerm.constants, outTensors) )

  # Normal ordering for sfExOps
  if True:

    # Make separate lists of the spin free excitation operators and other tensors
    #sfExOp_list = []
    #other_list  = []
    Op_list  = []
    for t in inTerm.tensors:
      Op_list.append(t.copy())
#      if isinstance(t, sfExOp):
#        sfExOp_list.append(t.copy())
#      else:
#        other_list.append(t.copy())

#    # Initialize n, the number of remaining spin free excitation operators
#    n = len(sfExOp_list)
#    if n != 1:
#        raise RuntimeError, "terms should have single sfExOp"

    # 
    iter_output_terms  = []
    ncore = 0
    nact  = 0
    nvirt = 0
    idx_dict = {}
    for t in Op_list:
        for idx in t.indices: 
            if idx.name not in idx_dict:
                #for typ in idx.indType:
                if len(idx.indType)==1:
                    typ = idx.indType[0]
                    if typ   == options.core_type:
                        idx_dict[idx.name] = 'i%d'%(ncore)
                        ncore += 1 
                    elif typ == options.active_type:
                        idx_dict[idx.name] = 'r%d'%(nact)
                        nact += 1 
                    elif typ == options.virtual_type:
                        idx_dict[idx.name] = 'a%d'%(nvirt)
                        nvirt += 1 

    for t in Op_list:
        idx_t_dict = {}
        ntot  = 0
        for idx in t.indices:
            if idx.name in idx_dict:
                idx.name = idx_dict[idx.name]
            elif idx.name in idx_t_dict:
                idx.name = idx_t_dict[idx.name]
            else:
                idx_t_dict[idx.name] = 'm%d'%(ntot)
                idx.name = idx_t_dict[idx.name]
                ntot += 1 

    outTerm = term(inTerm.numConstant, inTerm.constants, [Op_list[2]]+Op_list[0:2]+Op_list[3:5])

  else:
    raise RuntimeError, "HFFermi failed to choose what to do."

  return outTerm


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

def weight_factor(inTerm):
  """
  2. consider weight of spatial diagram   *** need to implement symmetric weight factor (1/2)
  """ 
  

#  if options.verbose:
#    print "converting to normal order:  %s" %(str(inTerm))

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError, "inTerm must be of class term"

#  # determine what types of operators the term contains
#  has_creDesOps = False
#  has_sfExOps = False
#  for t in inTerm.tensors:
#    if isinstance(t, creOp) or isinstance(t, desOp):
#      has_creDesOps = True
#    elif isinstance(t, sfExOp):
#      has_sfExOps = True
#
#  # If not spin free excitation operators,
#  # raise an error
#  if not has_sfExOps:
#    raise RuntimeError, "Present version of HFFermi can treat only sfExOp"
#
#  if has_creDesOps and has_sfExOps:
#    raise RuntimeError, "Present version of HFFermi can treat only sfExOp"

  # if the term is already normal ordered, return it unchanged
  elif not inTerm.isNormalOrdered():
    raise RuntimeError, "HFFermi requires normal ordered operators"

#  # Normal ordering for creOp/desOp
#  elif has_creDesOps:
#
#    # Separate the cre/des operators from other tensors
#    ops = []
#    nonOps = []
#    for t in inTerm.tensors:
#      if isinstance(t, creOp) or isinstance(t, desOp):
#        ops.append(t.copy())
#      else:
#        nonOps.append(t.copy())
#
#    # Generate all contraction pairs
#    contractionPairs = []
#    for i in range(len(ops)):
#      iTerm = ops[i]
#      for j in range(i+1,len(ops)):
#        jTerm = ops[j]
#        if isinstance(iTerm, desOp) and isinstance(jTerm, creOp):
#          contractionPairs.append((i,j));
#    #print "contractionPairs\n", contractionPairs
#
#    # Determine maximum contraction order
#    creCount = 0
#    maxConOrder = 0
#    for i in range(len(ops)-1,-1,-1):
#      iTerm = ops[i]
#      if isinstance(iTerm, creOp):
#        creCount +=1
#      elif isinstance(iTerm, desOp) and creCount > 0:
#        maxConOrder += 1
#        creCount -= 1
#    del(creCount,iTerm)
#
#    # Generate all contractions
#    contractions = []
#    for i in range(maxConOrder+1):
#      subCons = makeTuples(i,contractionPairs)
#      j = 0
#      while j < len(subCons):
#        creOpTags = []
#        desOpTags = []
#        for k in range(i):
#          creOpTags.append(subCons[j][k][1])
#          desOpTags.append(subCons[j][k][0])
#        if allDifferent(creOpTags) and allDifferent(desOpTags):
#          j += 1
#        else:
#          del(subCons[j])
#      for j in range(len(subCons)):
#        contractions.append(subCons[j])
#    del(subCons,creOpTags,desOpTags,contractionPairs)
#    #print "contractions:\n", contractions
#
#    # For each contraction, generate the resulting term
#    outTerms = []
#    for contraction in contractions:
#      conSign = 1
#      deltaFuncs = []
#      conIndeces = []
#      subOpString = []
#      subOpString.extend(ops)
#      for conPair in contraction:
#        index1 = ops[conPair[0]].indices[0]
#        index2 = ops[conPair[1]].indices[0]
#        deltaFuncs.append(kroneckerDelta([index1,index2]))
#        subOpString[conPair[0]] = 'contracted'
#        subOpString[conPair[1]] = 'contracted'
#        for q in subOpString[conPair[0]+1:conPair[1]]:
#          if not (q is 'contracted'):
#            conSign *= -1
#      i = 0
#      while i < len(subOpString):
#        if subOpString[i] is 'contracted':
#          del(subOpString[i])
#        else:
#          i += 1
#      (sortSign,sortedOps) = sortOps(subOpString)
#      totalSign = conSign * sortSign
#      outTensors = []
#      outTensors.extend(nonOps)
#      outTensors.extend(deltaFuncs)
#      outTensors.extend(sortedOps)
#      outTerms.append( term(totalSign * inTerm.numConstant, inTerm.constants, outTensors) )

  # Normal ordering for sfExOps
  if True:

    # Make separate lists of the spin free excitation operators and other tensors
    #sfExOp_list = []
    #other_list  = []
    Op_list  = []
    for t in inTerm.tensors:
      Op_list.append(t.copy())
#      if isinstance(t, sfExOp):
#        sfExOp_list.append(t.copy())
#      else:
#        other_list.append(t.copy())

#    # Initialize n, the number of remaining spin free excitation operators
#    n = len(sfExOp_list)
#    if n != 1:
#        raise RuntimeError, "terms should have single sfExOp"

    # 
    nloops = 0
    for i in range(len(Op_list)-1):
        ti = Op_list[i]
        orderi = len(ti.indices)/2
        for io in range(orderi):
            idx1i = ti.indices[io] 
            idx2i = ti.indices[io+orderi] 

            for j in range(len(Op_list)-i-1):
                tj = Op_list[j+i+1]
                orderj = len(tj.indices)/2
                for jo in range(orderj):
                    idx1j = tj.indices[jo] 
                    idx2j = tj.indices[jo+orderj] 
                    if idx1i == idx1j and idx2i == idx2j:
                        nloops += 1

    scal = 2**nloops
    outTerm = term(scal*inTerm.numConstant, inTerm.constants, Op_list)

  else:
    raise RuntimeError, "HFFermi failed to choose what to do."

  return outTerm



