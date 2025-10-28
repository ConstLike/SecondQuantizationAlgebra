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

def vacuumFermi(inTerm):
  ""

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError("inTerm must be of class term")

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
  if has_sfExOps:
    raise RuntimeError("Present version of vacuumFermi can treat only CreDesOps")

  # if the term is already normal ordered, return it unchanged
  elif not inTerm.isNormalOrdered():
    raise RuntimeError("vacuumFermi requires normal ordered operators")

  # Normal ordering for sfExOps
  elif has_creDesOps:

    const = inTerm.numConstant
    outTerm = term(0.0, [], inTerm.tensors)

  else:
    outTerm = inTerm

  return outTerm

def rmOverlap(inTerm):
  ""

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError("inTerm must be of class term")

  # determine what types of operators the term contains
  has_kronecker = False
  for t in inTerm.tensors:
    if isinstance(t, kroneckerDelta):
      has_kronecker = True 

  # Normal ordering for sfExOps
  if has_kronecker :
    outTerm = term(0.0, [], inTerm.tensors)
  else:
    outTerm = inTerm

  return outTerm

def Pauli(inTerm):
  ""

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError("inTerm must be of class term")

  # determine what types of operators the term contains
  has_creDesOps = False
  has_sfExOps   = False
  has_kronecker = False
  for t in inTerm.tensors:
    if isinstance(t, creOp) or isinstance(t, desOp):
      has_creDesOps = True
    elif isinstance(t, sfExOp):
      has_sfExOps   = True
    elif isinstance(t, kroneckerDelta):
      has_kronecker = True 


#  has_XG_XC2V = False
#  for t in inTerm.tensors:
#    if t.name == 'X_G' or t.name == 'X*_G' or t.name == 'X_C2V1' or t.name == 'X*_C2V1' or  t.name == 'X_C1V1' or t.name == 'X*_C1V1' or t.name == 'X_C2V2' or t.name == 'X*_C2V2':
#      has_XG_XC2V = True 

#  has_t2Wk = False
#  for t in inTerm.tensors:
#    if t.name == 't2Wk':
#      has_t2Wk = True

  # If not spin free excitation operators,
  # raise an error
  #if has_sfExOps or has_creDesOps or has_kronecker:
  #  raise RuntimeError, "Present version of Pauli only after contract"

  # if the term is already normal ordered, return it unchanged
  if not inTerm.isNormalOrdered():
    raise RuntimeError("HFFermi requires normal ordered operators")

  # Make separate lists of the spin free excitation operators and other tensors
  X_list  = []
  for t in inTerm.tensors:
    if t.name[0] == 'X':
      X_list.append(t.copy())

  const = inTerm.numConstant

  for x in X_list:
      if x.indices[0].name == x.indices[2].name: const = 0.0 
      if x.indices[1].name == x.indices[3].name: const = 0.0 

  outTerm = term(const, [], inTerm.tensors)

#  # Normal ordering for sfExOps
#  elif has_t2Wk:
#
#    # Make separate lists of the spin free excitation operators and other tensors
#    sfExOp_list = []
#    kroneckerDelta_list  = []
#    t2Wk_list  = []
#    t2Wb_list  = []
#    for t in inTerm.tensors:
#      if isinstance(t, sfExOp):
#        sfExOp_list.append(t.copy())
#      elif isinstance(t, kroneckerDelta):
#        kroneckerDelta_list.append(t.copy())
#      elif t.name == 't2Wk':
#        t2Wk_list.append(t.copy())
#      elif t.name == 't2Wb':
#        t2Wb_list.append(t.copy())
#      else:
#        raise RuntimeError, "terms should be sfExOp, kroneckerDelta, t2Wk, or t2Wb"
#
#    # Initialize n, the number of remaining spin free excitation operators
#    n = len(sfExOp_list)
#    if n != 1:
#        raise RuntimeError, "terms should have single sfExOp"
#
#    has_active = False
#    for idx in sfExOp_list[0].indices: 
#        for typ in idx.indType:
#            if typ == options.active_type:
#                has_active = True 
#
#    if not has_active:
#        const = inTerm.numConstant
#        order = sfExOp_list[0].order
#        for idx in sfExOp_list[0].indices: 
#            for typ in idx.indType:
#                #print typ[0]
#                if typ == options.virtual_type:
#                    const = 0.0 
#    
#        deltaFuncs = []
#        sfExOp_tmp= sfExOp_list[0]
#        order = sfExOp_tmp.order
#        if const == inTerm.numConstant:
#            a0_flag = False 
#            a2_flag = False
#            a1_flag = False
#    
#            # a0 
#            for d in kroneckerDelta_list:
#                if d.indices[0].name == "a0":
#                    a_delta = d.indices[1].name
#                elif d.indices[1].name == "a0":
#                    a_delta = d.indices[0].name
#    
#            for wb in t2Wb_list:
#                for i in range(3):
#                    if wb.indices[i+3].name == a_delta:
#                        i_t2Wb = wb.indices[i].name
#    
#            E_spin = sfExOp_tmp.spin
#            for i in range(order):
#                if sfExOp_tmp.indices[i].name   == "i0":
#                    spin_1 = E_spin[0][i]
#                if sfExOp_tmp.indices[i+3].name == i_t2Wb:
#                    spin_2 = E_spin[1][i]
#    
#            if spin_1 == spin_2: a0_flag = True
#    
#            # a2 
#            for d in kroneckerDelta_list:
#                if d.indices[0].name == "a2":
#                    a_delta = d.indices[1].name
#                elif d.indices[1].name == "a2":
#                    a_delta = d.indices[0].name
#    
#            for wb in t2Wb_list:
#                for i in range(3):
#                    if wb.indices[i+3].name == a_delta:
#                        i_t2Wb = wb.indices[i].name
#    
#            E_spin = sfExOp_tmp.spin
#            for i in range(order):
#                if sfExOp_tmp.indices[i].name   == "i2":
#                    spin_1 = E_spin[0][i]
#                if sfExOp_tmp.indices[i+3].name == i_t2Wb:
#                    spin_2 = E_spin[1][i]
#    
#            if spin_1 == spin_2: a2_flag = True
#    
#            # a1 
#            for d in kroneckerDelta_list:
#                if d.indices[0].name == "a1":
#                    a_delta = d.indices[1].name
#                elif d.indices[1].name == "a1":
#                    a_delta = d.indices[0].name
#    
#            for wb in t2Wb_list:
#                for i in range(3):
#                    if wb.indices[i+3].name == a_delta:
#                        i_t2Wb = wb.indices[i].name
#    
#            E_spin = sfExOp_tmp.spin
#            for i in range(order):
#                if sfExOp_tmp.indices[i].name   == "i1":
#                    spin_1 = E_spin[0][i]
#                if sfExOp_tmp.indices[i+3].name == i_t2Wb:
#                    spin_2 = E_spin[1][i]
#    
#            if spin_1 == spin_2: a1_flag = True
#    
#            if a1_flag is False:
#                if a0_flag is True or a2_flag is True: const *= 2
#            else:
#                if a0_flag is True and a2_flag is True: const *= 4 
#                if a0_flag is False and a2_flag is False: const *= 2 
#    
#    
#            for i in range(order):
#                index1 = sfExOp_list[0].indices[i]
#                index2 = sfExOp_list[0].indices[i+order]
#                deltaFuncs.append(kroneckerDelta([index1,index2])) 
#            #outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + t_list + W_list)
#            outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + t2Wk_list + t2Wb_list)
#    
#        else:
#            outTerm = term(const, inTerm.constants, deltaFuncs + t2Wk_list + t2Wb_list + sfExOp_list)
#            #outTerm = term(const, inTerm.constants, deltaFuncs + [t_list[0]] + W_list + [t_list[1]] + sfExOp_list)
#            #outTerm = term(const, inTerm.constants, other_list + sfExOp_list)
#
#    elif has_active:
#
#        const = inTerm.numConstant
#        order = sfExOp_list[0].order
#        for idx in sfExOp_list[0].indices: 
#            for typ in idx.indType:
#                #print typ[0]
#                if typ == options.virtual_type:
#                    const = 0.0 
#
#        deltaFuncs = []
#        sfExOp_tmp= sfExOp_list[0].copy()
#        order_tmp = sfExOp_list[0].order
#        E_spin = sfExOp_list[0].spin
#        E_index= sfExOp_list[0].indices
#        if const == inTerm.numConstant:
#            Dup = True 
#            while Dup:
#                Dup = False 
#                for i in range(order_tmp): 
#                    idxi = sfExOp_tmp.indices[order_tmp-i-1] 
#                    name_idxi = idxi.name
#                    type_idxi = idxi.indType
#                    for j in range(order_tmp): 
#                        idxj = sfExOp_tmp.indices[j+order_tmp] 
#                        name_idxj = idxj.name
#                        type_idxj = idxj.indType
#                        if type_idxi[0] == options.core_type and type_idxj[0] == options.core_type:
#                            Dup = True
#                            deltaFuncs.append(kroneckerDelta([idxi,idxj]))
#                            indices_tmp = []
#                            spin_tmp = [[],[]]
#                            for k in range(len(sfExOp_tmp.indices)): 
#                                if k != order_tmp-i-1:
#                                    if k != j+order_tmp:
#                                        indices_tmp.append(sfExOp_tmp.indices[k].copy())
#                                        if k <  order_tmp: spin_tmp[0].append(sfExOp_tmp.spin[0][k])
#                                        if k >= order_tmp: spin_tmp[1].append(sfExOp_tmp.spin[1][k-order_tmp])
#                            sfExOp_tmp = sfExOp(indices_tmp, spin=spin_tmp) 
#                            const *= (-1)**(j+order_tmp-(order_tmp-i-1)+1) 
#                            order_tmp = sfExOp_tmp.order
#                            #print i, j, name_idxi, name_idxj, order_tmp, sfExOp_tmp, sfExOp_tmp.spin
#                            break
#                    if Dup: break
#
#            # apply active spin
#            print sfExOp_tmp.spin
#            print E_spin
#            for i in range(len(sfExOp_tmp.spin[0])):
#                [sfExOp_tmp.spin[1][i] if x==sfExOp_tmp.spin[0][i] else x for x in E_spin[0]] 
#                [sfExOp_tmp.spin[1][i] if x==sfExOp_tmp.spin[0][i] else x for x in E_spin[1]] 
#            print E_spin
#
#            sfExOp_core_tmp= sfExOp_list[0].copy()
#            order_tmp      = sfExOp_list[0].order
#            indices_tmp    = []
#            spin_tmp       = [[],[]]
#            for k, idx in enumerate(sfExOp_core_tmp.indices):
#                if idx.indType[0] == options.core_type:
#                    indices_tmp.append(idx) 
#                    if k <  order_tmp: spin_tmp[0].append(sfExOp_core_tmp.spin[0][k])
#                    if k >= order_tmp: spin_tmp[1].append(sfExOp_core_tmp.spin[1][k-order_tmp])
#
#            sfExOp_core_tmp = sfExOp(indices_tmp, spin=spin_tmp) 
#            order_tmp = sfExOp_core_tmp.order
#
#            # scaling    
#            a0_flag = False 
#            a2_flag = False
#            a1_flag = False
#    
#            # a0 
#            if t2Wk_list[0].indices[0].indType[0] == options.core_type:
#                for d in kroneckerDelta_list:
#                    if d.indices[0].name == "a0":
#                        a_delta = d.indices[1].name
#                    elif d.indices[1].name == "a0":
#                        a_delta = d.indices[0].name
#        
#                for wb in t2Wb_list:
#                    for i in range(3):
#                        if wb.indices[i+3].name == a_delta:
#                            i_t2Wb = wb.indices[i].name
#        
#                for d in deltaFuncs:
#                    if sfExOp_core_tmp.indices[i].name   == "i0":
#                        spin_1 = E_spin[0][i]
#    
#                for i in range(order_tmp):
#                    if sfExOp_core_tmp.indices[i].name   == "i0":
#                        spin_1 = E_spin[0][i]
#                    if sfExOp_core_tmp.indices[i+3].name == i_t2Wb:
#                        spin_2 = E_spin[1][i]
#        
#                if spin_1 == spin_2: a0_flag = True
#    
#            # a2 
#            if t2Wk_list[0].indices[2].indType[0] == options.core_type:
#                for d in kroneckerDelta_list:
#                    if d.indices[0].name == "a2":
#                        a_delta = d.indices[1].name
#                    elif d.indices[1].name == "a2":
#                        a_delta = d.indices[0].name
#        
#                for wb in t2Wb_list:
#                    for i in range(3):
#                        if wb.indices[i+3].name == a_delta:
#                            i_t2Wb = wb.indices[i].name
#        
#                E_spin = sfExOp_core_tmp.spin
#                for i in range(order_tmp):
#                    if sfExOp_core_tmp.indices[i].name   == "i2":
#                        spin_1 = E_spin[0][i]
#                    if sfExOp_core_tmp.indices[i+3].name == i_t2Wb:
#                        spin_2 = E_spin[1][i]
#        
#                if spin_1 == spin_2: a2_flag = True
#    
#            # a1 
#            if t2Wk_list[0].indices[1].indType[0] == options.core_type:
#                for d in kroneckerDelta_list:
#                    if d.indices[0].name == "a1":
#                        a_delta = d.indices[1].name
#                    elif d.indices[1].name == "a1":
#                        a_delta = d.indices[0].name
#        
#                for wb in t2Wb_list:
#                    for i in range(3):
#                        if wb.indices[i+3].name == a_delta:
#                            i_t2Wb = wb.indices[i].name
#        
#                E_spin = sfExOp_core_tmp.spin
#                for i in range(order_tmp):
#                    if sfExOp_core_tmp.indices[i].name   == "i1":
#                        spin_1 = E_spin[0][i]
#                    if sfExOp_core_tmp.indices[i+3].name == i_t2Wb:
#                        spin_2 = E_spin[1][i]
#        
#                if spin_1 == spin_2: a1_flag = True
#    
#            if a1_flag is False:
#                if a0_flag is True  or  a2_flag is True:  const *= 2
#            else:
#                if a0_flag is True  and a2_flag is True:  const *= 4 
#                if a0_flag is False and a2_flag is False: const *= 2 
#    
##    
##            for i in range(order_tmp):
##                index1 = sfExOp_list[0].indices[i]
##                index2 = sfExOp_list[0].indices[i+order_tmp]
##                deltaFuncs.append(kroneckerDelta([index1,index2])) 
#            #outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + t_list + W_list)
#            outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + t2Wk_list + t2Wb_list + [sfExOp_tmp])
#    
#        else:
#            outTerm = term(const, inTerm.constants, deltaFuncs + t2Wk_list + t2Wb_list + sfExOp_list)
#            #outTerm = term(const, inTerm.constants, deltaFuncs + [t_list[0]] + W_list + [t_list[1]] + sfExOp_list)
#            #outTerm = term(const, inTerm.constants, other_list + sfExOp_list)

#  else:
#    raise RuntimeError, "Pauli failed to choose what to do."

  return outTerm

def rmTags(inTerm):
  ""

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError("inTerm must be of class term")

  # determine what types of operators the term contains
  has_creDesOps = False
  has_sfExOps   = False
  has_kronecker = False
  for t in inTerm.tensors:
    if isinstance(t, creOp) or isinstance(t, desOp):
      has_creDesOps = True
    elif isinstance(t, sfExOp):
      has_sfExOps   = True
    elif isinstance(t, kroneckerDelta):
      has_kronecker = True 

#  if has_sfExOps or has_creDesOps or has_kronecker:
#    raise RuntimeError, "Present version of Pauli only after contract"
  if not inTerm.isNormalOrdered():
    raise RuntimeError("HFFermi requires normal ordered operators")

  new_tensors = []  
  for t in inTerm.tensors:
    new_indices = []
    for ind in t.indices:
        indn = index(ind.name, [], ind.isSummed)
        new_indices.append(indn)

    tn = tensor(t.name, new_indices, t.symmetries) 
    new_tensors.append(tn)

  outTerm = term(inTerm.numConstant, [], new_tensors)

  return outTerm


def HFFermi(inTerm):
  ""

#  if options.verbose:
#    print "converting to normal order:  %s" %(str(inTerm))

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError("inTerm must be of class term")

  # determine what types of operators the term contains
  has_creDesOps = False
  has_sfExOps = False
  for t in inTerm.tensors:
    if isinstance(t, creOp) or isinstance(t, desOp):
      has_creDesOps = True
    elif isinstance(t, sfExOp):
      has_sfExOps = True


  has_W = False
  for t in inTerm.tensors:
    if t.name == 'W':
      has_W = True 

  has_t2Wk = False
  for t in inTerm.tensors:
    if t.name == 't2Wk':
      has_t2Wk = True

  # If not spin free excitation operators,
  # raise an error
  if not has_sfExOps:
    raise RuntimeError("Present version of HFFermi can treat only sfExOp")

  if has_creDesOps and has_sfExOps:
    raise RuntimeError("Present version of HFFermi can treat only sfExOp")

  # if the term is already normal ordered, return it unchanged
  elif not inTerm.isNormalOrdered():
    raise RuntimeError("HFFermi requires normal ordered operators")


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
  elif has_W:

    # Make separate lists of the spin free excitation operators and other tensors
    sfExOp_list = []
    kroneckerDelta_list  = []
    W_list  = []
    t_list  = []
    for t in inTerm.tensors:
      if isinstance(t, sfExOp):
        sfExOp_list.append(t.copy())
      elif isinstance(t, kroneckerDelta):
        kroneckerDelta_list.append(t.copy())
      elif t.name == 'W':
        W_list.append(t.copy())
      elif t.name == 't2':
        t_list.append(t.copy())
      else:
        raise RuntimeError("terms should be sfExOp, kroneckerDelta, W, or t2")

    # Initialize n, the number of remaining spin free excitation operators
    n = len(sfExOp_list)
    if n != 1:
        raise RuntimeError("terms should have single sfExOp")

    const = inTerm.numConstant
    order = sfExOp_list[0].order
    for idx in sfExOp_list[0].indices: 
        for typ in idx.indType:
            #print typ[0]
    #TODO: when a,b exist in E, delta term should be remained
    #TODO: this version just ignore the term 
            if typ == options.virtual_type:
                const = 0.0 

    deltaFuncs = []
    sfExOp_tmp= sfExOp_list[0]
    order = sfExOp_tmp.order
    if const == inTerm.numConstant:
        # scaling
        # assign t2 spin from E spin
        # sym1 : list of indices of sigma, sym2 : list of indices of tau
        sym_list = []        
        for t in t_list:
            for i in range(2):   # 2 = len(t.indices)//2
                sym_list.append( [t.indices[i].name, t.indices[i+2].name] )
        #print sym_list 

        for d in kroneckerDelta_list:
            sym_list.append( [d.indices[0].name, d.indices[1].name] )
        #print sym_list 

        for i in range(order):
            s1 = sfExOp_tmp.spin[0][i]
            j  = sfExOp_tmp.spin[1].index(s1)
            sym_list.append( [sfExOp_tmp.indices[i].name, sfExOp_tmp.indices[j].name] )
        #print sym_list 
            
        for w in W_list:
            for i in range(2):   # 2 = len(t.indices)//2
                sym_list.append( [w.indices[i].name, w.indices[i+2].name] )
        #print sym_list 

        def common_member(a, b):
            a_set = set(a)
            b_set = set(b)
            if (a_set & b_set):
                return True 
            else:
                return False 

        def common_list(sym_list, idx_c):
            no_comm = True 
            idx_nc = []
            for sym in sym_list:
                if common_member( idx_c, sym ):
                    idx_c  += sym
                    no_comm = False 
                else:
                    idx_nc += sym
            return idx_nc, no_comm

        for i, w in enumerate(W_list):
            w_s = [w.indices[0].name, w.indices[2].name]
            w_t = [w.indices[1].name, w.indices[3].name]
            idx_c  = w_s[:] 
            idx_nc = sym_list[:]
            no_comm = True 

            while no_comm:
                idx_nc, no_comm = common_list(idx_nc, idx_c)

            if not common_member ( w_t, idx_c ): const *= 2


        for i in range(order):
            index1 = sfExOp_list[0].indices[i]
            index2 = sfExOp_list[0].indices[i+order]
            deltaFuncs.append(kroneckerDelta([index1,index2])) 
        #outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + t_list + W_list)
        outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + [t_list[0]] + W_list + [t_list[1]])

    else:
        outTerm = term(const, inTerm.constants, deltaFuncs + t_list + W_list + sfExOp_list)
        #outTerm = term(const, inTerm.constants, deltaFuncs + [t_list[0]] + W_list + [t_list[1]] + sfExOp_list)
        #outTerm = term(const, inTerm.constants, other_list + sfExOp_list)

  # Normal ordering for sfExOps
  elif has_t2Wk:

    # Make separate lists of the spin free excitation operators and other tensors
    sfExOp_list = []
    kroneckerDelta_list  = []
    t2Wk_list  = []
    t2Wb_list  = []
    for t in inTerm.tensors:
      if isinstance(t, sfExOp):
        sfExOp_list.append(t.copy())
      elif isinstance(t, kroneckerDelta):
        kroneckerDelta_list.append(t.copy())
      elif t.name == 't2Wk':
        t2Wk_list.append(t.copy())
      elif t.name == 't2Wb':
        t2Wb_list.append(t.copy())
      else:
        raise RuntimeError("terms should be sfExOp, kroneckerDelta, t2Wk, or t2Wb")

    # Initialize n, the number of remaining spin free excitation operators
    n = len(sfExOp_list)
    if n != 1:
        raise RuntimeError("terms should have single sfExOp")

    has_active = False
    for idx in sfExOp_list[0].indices: 
        for typ in idx.indType:
            if typ == options.active_type:
                has_active = True 

    if not has_active:
        const = inTerm.numConstant
        order = sfExOp_list[0].order
        for idx in sfExOp_list[0].indices: 
            for typ in idx.indType:
                #print typ[0]
                if typ == options.virtual_type:
                    const = 0.0 
    
        deltaFuncs = []
        sfExOp_tmp= sfExOp_list[0]
        order = sfExOp_tmp.order
        if const == inTerm.numConstant:
            a0_flag = False 
            a2_flag = False
            a1_flag = False
    
            # a0 
            for d in kroneckerDelta_list:
                if d.indices[0].name == "a0":
                    a_delta = d.indices[1].name
                elif d.indices[1].name == "a0":
                    a_delta = d.indices[0].name
    
            for wb in t2Wb_list:
                for i in range(3):
                    if wb.indices[i+3].name == a_delta:
                        i_t2Wb = wb.indices[i].name
    
            E_spin = sfExOp_tmp.spin
            for i in range(order):
                if sfExOp_tmp.indices[i].name   == "i0":
                    spin_1 = E_spin[0][i]
                if sfExOp_tmp.indices[i+3].name == i_t2Wb:
                    spin_2 = E_spin[1][i]
    
            if spin_1 == spin_2: a0_flag = True
    
            # a2 
            for d in kroneckerDelta_list:
                if d.indices[0].name == "a2":
                    a_delta = d.indices[1].name
                elif d.indices[1].name == "a2":
                    a_delta = d.indices[0].name
    
            for wb in t2Wb_list:
                for i in range(3):
                    if wb.indices[i+3].name == a_delta:
                        i_t2Wb = wb.indices[i].name
    
            E_spin = sfExOp_tmp.spin
            for i in range(order):
                if sfExOp_tmp.indices[i].name   == "i2":
                    spin_1 = E_spin[0][i]
                if sfExOp_tmp.indices[i+3].name == i_t2Wb:
                    spin_2 = E_spin[1][i]
    
            if spin_1 == spin_2: a2_flag = True
    
            # a1 
            for d in kroneckerDelta_list:
                if d.indices[0].name == "a1":
                    a_delta = d.indices[1].name
                elif d.indices[1].name == "a1":
                    a_delta = d.indices[0].name
    
            for wb in t2Wb_list:
                for i in range(3):
                    if wb.indices[i+3].name == a_delta:
                        i_t2Wb = wb.indices[i].name
    
            E_spin = sfExOp_tmp.spin
            for i in range(order):
                if sfExOp_tmp.indices[i].name   == "i1":
                    spin_1 = E_spin[0][i]
                if sfExOp_tmp.indices[i+3].name == i_t2Wb:
                    spin_2 = E_spin[1][i]
    
            if spin_1 == spin_2: a1_flag = True
    
            if a1_flag is False:
                if a0_flag is True or a2_flag is True: const *= 2
            else:
                if a0_flag is True and a2_flag is True: const *= 4 
                if a0_flag is False and a2_flag is False: const *= 2 
    
    
            for i in range(order):
                index1 = sfExOp_list[0].indices[i]
                index2 = sfExOp_list[0].indices[i+order]
                deltaFuncs.append(kroneckerDelta([index1,index2])) 
            #outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + t_list + W_list)
            outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + t2Wk_list + t2Wb_list)
    
        else:
            outTerm = term(const, inTerm.constants, deltaFuncs + t2Wk_list + t2Wb_list + sfExOp_list)
            #outTerm = term(const, inTerm.constants, deltaFuncs + [t_list[0]] + W_list + [t_list[1]] + sfExOp_list)
            #outTerm = term(const, inTerm.constants, other_list + sfExOp_list)

    elif has_active:

        const = inTerm.numConstant
        order = sfExOp_list[0].order
        for idx in sfExOp_list[0].indices: 
            for typ in idx.indType:
                #print typ[0]
                if typ == options.virtual_type:
                    const = 0.0 

        deltaFuncs = []
        sfExOp_tmp= sfExOp_list[0].copy()
        order_tmp = sfExOp_list[0].order
        E_spin = sfExOp_list[0].spin
        E_index= sfExOp_list[0].indices
        if const == inTerm.numConstant:
            Dup = True 
            while Dup:
                Dup = False 
                for i in range(order_tmp): 
                    idxi = sfExOp_tmp.indices[order_tmp-i-1] 
                    name_idxi = idxi.name
                    type_idxi = idxi.indType
                    for j in range(order_tmp): 
                        idxj = sfExOp_tmp.indices[j+order_tmp] 
                        name_idxj = idxj.name
                        type_idxj = idxj.indType
                        if type_idxi[0] == options.core_type and type_idxj[0] == options.core_type:
                            Dup = True
                            deltaFuncs.append(kroneckerDelta([idxi,idxj]))
                            indices_tmp = []
                            spin_tmp = [[],[]]
                            for k in range(len(sfExOp_tmp.indices)): 
                                if k != order_tmp-i-1:
                                    if k != j+order_tmp:
                                        indices_tmp.append(sfExOp_tmp.indices[k].copy())
                                        if k <  order_tmp: spin_tmp[0].append(sfExOp_tmp.spin[0][k])
                                        if k >= order_tmp: spin_tmp[1].append(sfExOp_tmp.spin[1][k-order_tmp])
                            sfExOp_tmp = sfExOp(indices_tmp, spin=spin_tmp) 
                            const *= (-1)**(j+order_tmp-(order_tmp-i-1)+1) 
                            order_tmp = sfExOp_tmp.order
                            #print i, j, name_idxi, name_idxj, order_tmp, sfExOp_tmp, sfExOp_tmp.spin
                            break
                    if Dup: break

            # apply active spin
            print(sfExOp_tmp.spin)
            print(E_spin)
            for i in range(len(sfExOp_tmp.spin[0])):
                [sfExOp_tmp.spin[1][i] if x==sfExOp_tmp.spin[0][i] else x for x in E_spin[0]] 
                [sfExOp_tmp.spin[1][i] if x==sfExOp_tmp.spin[0][i] else x for x in E_spin[1]] 
            print(E_spin)

            sfExOp_core_tmp= sfExOp_list[0].copy()
            order_tmp      = sfExOp_list[0].order
            indices_tmp    = []
            spin_tmp       = [[],[]]
            for k, idx in enumerate(sfExOp_core_tmp.indices):
                if idx.indType[0] == options.core_type:
                    indices_tmp.append(idx) 
                    if k <  order_tmp: spin_tmp[0].append(sfExOp_core_tmp.spin[0][k])
                    if k >= order_tmp: spin_tmp[1].append(sfExOp_core_tmp.spin[1][k-order_tmp])

            sfExOp_core_tmp = sfExOp(indices_tmp, spin=spin_tmp) 
            order_tmp = sfExOp_core_tmp.order

            # scaling    
            a0_flag = False 
            a2_flag = False
            a1_flag = False
    
            # a0 
            if t2Wk_list[0].indices[0].indType[0] == options.core_type:
                for d in kroneckerDelta_list:
                    if d.indices[0].name == "a0":
                        a_delta = d.indices[1].name
                    elif d.indices[1].name == "a0":
                        a_delta = d.indices[0].name
        
                for wb in t2Wb_list:
                    for i in range(3):
                        if wb.indices[i+3].name == a_delta:
                            i_t2Wb = wb.indices[i].name
        
                for d in deltaFuncs:
                    if sfExOp_core_tmp.indices[i].name   == "i0":
                        spin_1 = E_spin[0][i]
    
                for i in range(order_tmp):
                    if sfExOp_core_tmp.indices[i].name   == "i0":
                        spin_1 = E_spin[0][i]
                    if sfExOp_core_tmp.indices[i+3].name == i_t2Wb:
                        spin_2 = E_spin[1][i]
        
                if spin_1 == spin_2: a0_flag = True
    
            # a2 
            if t2Wk_list[0].indices[2].indType[0] == options.core_type:
                for d in kroneckerDelta_list:
                    if d.indices[0].name == "a2":
                        a_delta = d.indices[1].name
                    elif d.indices[1].name == "a2":
                        a_delta = d.indices[0].name
        
                for wb in t2Wb_list:
                    for i in range(3):
                        if wb.indices[i+3].name == a_delta:
                            i_t2Wb = wb.indices[i].name
        
                E_spin = sfExOp_core_tmp.spin
                for i in range(order_tmp):
                    if sfExOp_core_tmp.indices[i].name   == "i2":
                        spin_1 = E_spin[0][i]
                    if sfExOp_core_tmp.indices[i+3].name == i_t2Wb:
                        spin_2 = E_spin[1][i]
        
                if spin_1 == spin_2: a2_flag = True
    
            # a1 
            if t2Wk_list[0].indices[1].indType[0] == options.core_type:
                for d in kroneckerDelta_list:
                    if d.indices[0].name == "a1":
                        a_delta = d.indices[1].name
                    elif d.indices[1].name == "a1":
                        a_delta = d.indices[0].name
        
                for wb in t2Wb_list:
                    for i in range(3):
                        if wb.indices[i+3].name == a_delta:
                            i_t2Wb = wb.indices[i].name
        
                E_spin = sfExOp_core_tmp.spin
                for i in range(order_tmp):
                    if sfExOp_core_tmp.indices[i].name   == "i1":
                        spin_1 = E_spin[0][i]
                    if sfExOp_core_tmp.indices[i+3].name == i_t2Wb:
                        spin_2 = E_spin[1][i]
        
                if spin_1 == spin_2: a1_flag = True
    
            if a1_flag is False:
                if a0_flag is True  or  a2_flag is True:  const *= 2
            else:
                if a0_flag is True  and a2_flag is True:  const *= 4 
                if a0_flag is False and a2_flag is False: const *= 2 
    
#    
#            for i in range(order_tmp):
#                index1 = sfExOp_list[0].indices[i]
#                index2 = sfExOp_list[0].indices[i+order_tmp]
#                deltaFuncs.append(kroneckerDelta([index1,index2])) 
            #outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + t_list + W_list)
            outTerm = term(const, inTerm.constants, kroneckerDelta_list + deltaFuncs + t2Wk_list + t2Wb_list + [sfExOp_tmp])
    
        else:
            outTerm = term(const, inTerm.constants, deltaFuncs + t2Wk_list + t2Wb_list + sfExOp_list)
            #outTerm = term(const, inTerm.constants, deltaFuncs + [t_list[0]] + W_list + [t_list[1]] + sfExOp_list)
            #outTerm = term(const, inTerm.constants, other_list + sfExOp_list)

  else:
    raise RuntimeError("HFFermi failed to choose what to do.")

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
    raise TypeError("inTerm must be of class term")

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
    raise RuntimeError("HFFermi requires normal ordered operators")

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

#    # Op order determine
#    opidx = [0, 1, 2, 3]
#    for t in Op_list:
#        if t.name == 't2' and t.indices[0].indType == options.virtual_type

    #outTerm = term(inTerm.numConstant, inTerm.constants, [Op_list[2]]+Op_list[0:2]+Op_list[3:5])
    outTerm = term(inTerm.numConstant, inTerm.constants, Op_list)

  else:
    raise RuntimeError("HFFermi failed to choose what to do.")

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
    raise TypeError("inTerm must be of class term")

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
    raise RuntimeError("HFFermi requires normal ordered operators")

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
#                    print idx1i.name, idx2i.name, idx1j.name, idx2j.name 
                    if idx1i.name == idx1j.name and idx2i.name == idx2j.name or \
                       idx1i.name == idx2j.name and idx2i.name == idx1j.name:
                        nloops += 1



    scal = 2**nloops



    outTerm = term(scal*inTerm.numConstant, inTerm.constants, Op_list)

  else:
    raise RuntimeError("HFFermi failed to choose what to do.")

  return outTerm


def sort_noeqidx_terms(Terms_List):
  """
  """ 

  eq_l  = []
  neq_l = []
  for terms in Terms_List:
    neq = True
    idx_dict = {}
    for t in terms.tensors:
      for idx in t.indices:
        if idx.name not in idx_dict:
          idx_dict[idx.name]  = 1 
        else:
          idx_dict[idx.name] += 1 

    test = sum(1 for i in list(idx_dict.values()) if i >= 3) 

    if test == 0:
      neq_l.append(terms)
    else:          
      eq_l.append(terms)

  return ['eq'] + eq_l + ['noneq'] + neq_l 

def SpinFree(inTerm):
  ""
  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError("inTerm must be of class term")

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
    raise RuntimeError("Present version of HFFermi can treat only sfExOp")

  if has_creDesOps and has_sfExOps:
    raise RuntimeError("Present version of HFFermi can treat only sfExOp")

  # if the term is already normal ordered, return it unchanged
  elif not inTerm.isNormalOrdered():
    raise RuntimeError("HFFermi requires normal ordered operators")

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
        raise RuntimeError("terms should have single sfExOp")

    # obtain rm same index in Spin-Free Op 
    sfExOp_tmp = sfExOp_list[0]
    const_tmp = inTerm.numConstant
    order_tmp = sfExOp_list[0].order
    Dup = True 
    while Dup:
        Dup = False 
        for i in range(order_tmp): 
            idxi = sfExOp_tmp.indices[order_tmp-i-1] 
            name_idxi = idxi.name
            for j in range(order_tmp): 
                idxj = sfExOp_tmp.indices[j+order_tmp] 
                name_idxj = idxj.name
                if name_idxi == name_idxj:
                    Dup = True
                    indices_tmp = []
                    spin_tmp = [[],[]]
                    for k in range(len(sfExOp_tmp.indices)): 
                        if k != order_tmp-i-1:
                            if k != j+order_tmp:
                                indices_tmp.append(sfExOp_tmp.indices[k].copy())
                                if k <  order_tmp: spin_tmp[0].append(sfExOp_tmp.spin[0][k])
                                if k >= order_tmp: spin_tmp[1].append(sfExOp_tmp.spin[1][k])
                    sfExOp_tmp = sfExOp(indices_tmp, spin=spin_tmp) 
                    const_tmp *= (-1)**(j+order_tmp-(order_tmp-i-1)+1) 
                    order_tmp = sfExOp_tmp.order
                    #print i, j, name_idxi, name_idxj, order_tmp, sfExOp_tmp
                    break
            if Dup: break

    outTerm = term(const_tmp, inTerm.constants, other_list + [sfExOp_tmp])

  else:
    raise RuntimeError("HFFermi failed to choose what to do.")

  return outTerm


