import secondQuantizationAlgebra as sqa
import math 

sqa.options.verbose = False
# tags 
tag_alpha   = sqa.options.alpha_type
tag_beta    = sqa.options.beta_type

tag_core    = sqa.options.core_type
tag_active1 = sqa.options.active_type
tag_active2 = sqa.options.active_type2
tag_virtual = sqa.options.virtual_type
tag_total  = tag_core + tag_active1 + tag_active2 + tag_virtual

# indices for wavef 
dummy = True
ia = sqa.index('ia', [tag_core, tag_alpha],   dummy)
ib = sqa.index('ib', [tag_core, tag_beta ],   dummy)
o1a = sqa.index('o1a', [tag_active1, tag_alpha],dummy)
o1b = sqa.index('o1b', [tag_active1, tag_beta ],dummy)
o2a = sqa.index('o2a', [tag_active2, tag_alpha],dummy)
o2b = sqa.index('o2b', [tag_active2, tag_beta ],dummy)
aa = sqa.index('aa', [tag_virtual, tag_alpha],dummy)
ab = sqa.index('ab', [tag_virtual, tag_beta ],dummy)
# indices for Op 
dummy = False 
ta = sqa.index('ta', [tag_total, tag_alpha],  dummy)
tb = sqa.index('tb', [tag_total, tag_beta ],  dummy)
ua = sqa.index('ua', [tag_total, tag_alpha],  dummy)
ub = sqa.index('ub', [tag_total, tag_beta ],  dummy)

# 
ia_cre = sqa.creOp(ia)
ib_cre = sqa.creOp(ib)
o1a_cre = sqa.creOp(o1a)
o1b_cre = sqa.creOp(o1b)
o2a_cre = sqa.creOp(o2a)
o2b_cre = sqa.creOp(o2b)
aa_cre = sqa.creOp(aa)
ab_cre = sqa.creOp(ab)

ia_des = sqa.desOp(ia)
ib_des = sqa.desOp(ib)
o1a_des = sqa.desOp(o1a)
o1b_des = sqa.desOp(o1b)
o2a_des = sqa.desOp(o2a)
o2b_des = sqa.desOp(o2b)
aa_des = sqa.desOp(aa)
ab_des = sqa.desOp(ab)

ta_des = sqa.desOp(ta)
tb_des = sqa.desOp(tb)
ua_cre = sqa.creOp(ua)
ub_cre = sqa.creOp(ub)

t_ta_des = sqa.term(1.0,[],[ta_des])
t_tb_des = sqa.term(1.0,[],[tb_des])
t_ua_cre = sqa.term(1.0,[],[ua_cre])
t_ub_cre = sqa.term(1.0,[],[ub_cre])

Ms_p_ref_bra = [o2a_cre, o1a_cre, ib_cre, ia_cre]
Ms_m_ref_bra = [o2b_cre, o1b_cre, ib_cre, ia_cre]

Ms_p_ref_ket = [ia_des, ib_des, o1a_des, o2a_des]
Ms_m_ref_ket = [ia_des, ib_des, o1b_des, o2b_des]

def genCSF(string, S, Ms, bra=False, ket=False):
    assert len(string) == 4                      # MRSF C, O1, O2, V orbital space
    for occ in string:
        assert occ == 0 or occ == 1 or occ == 2  # a string of occupation number (0/1/2) in orbital
    assert bra or ket                            # select bra or ket
    assert not (bra and ket)                     # select bra or ket
    nc1 = string[0]
    nc2 = string[1]
    nv1 = string[2]
    nv2 = string[3]

    assert (S == 0 and Ms == 0) or (S == 1 and Ms == 0)
    sign = 1
    if S == 1 and Ms == 0: sign = -1

    #Ms=+1 ref Op= [ia_cre,ib_cre,ja_cre,aa_cre]
    #Ms=-1 ref Op= [ia_cre,ib_cre,jb_cre,ab_cre]
    #G
    if nc1 == 2 and nc2 == 2 and nv1 == 0 and nv2 == 0:
        CGcoeff = 1.0
        if ket:
            X = [sqa.tensor('X_G', [o1b,o2a], [])]
            Op= [o1b_cre, o1a_cre, ib_cre, ia_cre]
            terms = [sqa.term(CGcoeff,[],X+Op)]
        elif bra:
            X = [sqa.tensor('X*_G', [o1b,o2a], [])]
            Op= [ia_des, ib_des, o1a_des, o1b_des]
            terms = [sqa.term(CGcoeff,[],X+Op)]
    #L,R
    elif nc1 == 2 and nc2 == 1 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_LR', [o1b,o1a], [])]
            Op= [o2a_cre, o1b_cre, ib_cre, ia_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_LR', [o2b,o2a], [])]
            Op= [o2b_cre, o1a_cre, ib_cre, ia_cre]
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_LR', [o1b,o1a], [])]
            Op= [ia_des, ib_des, o1b_des, o2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_LR', [o2b,o2a], [])]
            Op= [ia_des, ib_des, o1a_des, o2b_des]
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))

    #CO1
    elif nc1 == 1 and nc2 == 2 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_CO1', [o1b,ia], [])]
            Op= [o2a_cre, o1a_cre, ib_cre, o1b_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_CO1', [o1b,ia], [])]
            Op= [o2b_cre, o1b_cre, o1a_cre, ia_cre]
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_CO1', [o1b,ia], [])]
            Op= [o1b_des, ib_des, o1a_des, o2a_des] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_CO1', [o1b,ia], [])]
            Op= [ia_des, o1a_des, o1b_des, o2b_des]
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
    #O2V
    elif nc1 == 2 and nc2 == 1 and nv1 == 0 and nv2 == 1:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_O2V', [ab,o2a], [])]
            Op= [ab_cre, o1a_cre, ib_cre, ia_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_O2V', [ab,o2a], [])]
            Op= [aa_cre, o1b_cre, ib_cre, ia_cre]
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O2V', [ab,o2a], [])]
            Op= [ia_des, ib_des, o1a_des, ab_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_O2V', [ab,o2a], [])]
            Op= [ia_des, ib_des, o1b_des, aa_des]
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))

#Ms_p_ref_ket = [o2a_cre, o1a_cre, ib_cre, ia_cre]
#Ms_m_ref_ket = [o2b_cre, o1b_cre, ib_cre, ia_cre]
#
#Ms_p_ref_bra = [ia_des, ib_des, o1a_des, o2a_des]
#Ms_m_ref_bra = [ia_des, ib_des, o1b_des, o2b_des]



    #O1O2
    elif nc1 == 2 and nc2 == 0 and nv1 == 2 and nv2 == 0:
        CGcoeff = 1.0
        if ket:
            X = [sqa.tensor('X_G', [o2b,o1a], [])]
            Op= [o2a_cre, o2b_cre, ib_cre, ia_cre] 
            terms = [sqa.term(CGcoeff,[],X+Op)]
        elif bra:
            X = [sqa.tensor('X*_G', [o2b,o1a], [])]
            Op= [ia_des, ib_des, o2b_des, o2a_des]
            terms = [sqa.term(CGcoeff,[],X+Op)]
    #CO2
    elif nc1 == 1 and nc2 == 1 and nv1 == 2 and nv2 == 0:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_CO2', [o2b,ia], [])]
            Op= [o2a_cre, o1a_cre, ib_cre, o2b_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_CO2', [o2b,ia], [])]
            Op= [o2b_cre, o1b_cre, o2a_cre, ia_cre] 
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_CO2', [o2b,ia], [])]
            Op= [o2b_des, ib_des, o1a_des, o2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_CO2', [o2b,ia], [])]
            Op= [ia_des, o2a_des, o1b_des, o2b_des]
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
    #O1V
    elif nc1 == 2 and nc2 == 0 and nv1 == 1 and nv2 == 1:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_O1V', [ab,o1a], [])]
            Op= [o2a_cre, ab_cre, ib_cre, ia_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_O1V', [ab,o1a], [])]
            Op= [o2b_cre, aa_cre, ib_cre, ia_cre] 
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O1V', [ab,o1a], [])]
            Op= [ia_des, ib_des, ab_des, o2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_O1V', [ab,o1a], [])]
            Op= [ia_des, ib_des, aa_des, o2b_des]
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
    #CV - spin contaminated one
    elif nc1 == 1 and nc2 == 1 and nv1 == 1 and nv2 == 1:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_CV', [ab,ia], [])]
            Op= [o2a_cre, o1a_cre, ib_cre, ab_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_CV', [ab,ia], [])]
            Op= [o2b_cre, o1b_cre, aa_cre, ia_cre] 
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_CV', [ab,ia], [])]
            Op= [ab_des, ib_des, o1a_des, o2a_des] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_CV', [ab,ia], [])]
            Op= [ia_des, aa_des, o1b_des, o2b_des]
            terms.append(sqa.term(sign*CGcoeff,[],X+Op))
    else:
        assert False 
    return terms 

def overlap(bra, ops, ket, dbg=False):
    print("bra:")
    for t in bra:
        print(t)
    print("ops:")
    for t in ops:
        print(t)
    print("ket:")
    for t in ket:
        print(t)
    print("<bra|ops|ket>:")

    if dbg: print("expand terms")
    ov_tmp = []
    if len(ops) != 0:
        for tb in bra:
            opt = tb
            for op in ops: 
                opt = sqa.multiplyTerms(opt, op)
            for tk in ket:
                ov_tmp.append(sqa.multiplyTerms(opt, tk))
    else: 
        for tb in bra:
            for tk in ket:
                ov_tmp.append(sqa.multiplyTerms(tb, tk))

    #    normal order form
    if dbg: print("normal Order for maximum contraction form")
    ov_tmp2= []
    i = 0
    for t in ov_tmp:
        i += 1
        #print "%d / %d "%(i, len(ov_tmp))
        term = sqa.normalOrder_maxcontraction( t )
        #for tn in term:
        #    print tn
        ov_tmp2 += term

    #TODO: remove terms by vacuumFermi in sqaSLee
    if dbg: print("remove delta_{tu} and delta")
    ov_tmp3= []
    for t in ov_tmp2:
        tn = sqa.vacuumFermi( t )
        #print tn 
        ov_tmp3.append( tn )
    if dbg: print("")

    if dbg: print("contract form")
    for t in ov_tmp3:
        t.contractDeltaFuncs()
        if dbg: print(t)
    if dbg: print("")
    
    if dbg: print("remove near zero terms")
    sqa.termChop(ov_tmp3)
    if dbg: 
        for t in ov_tmp3:
            print(t)
        print("")

#    if dbg: print "apply Pauli principle"
#    ov = []
#    for t in ov_tmp3:
#        tn = sqa.Pauli( t )
#        if dbg: print tn 
#        ov.append( tn )
#    if dbg:  print ""

    if dbg: print("rm overlap")
    ov2 = []
    for t in ov_tmp3:
        tn = sqa.rmOverlap( t )
        if dbg: print(tn) 
        ov2.append( tn )
    if dbg:  print("")

    if dbg: print("remove near zero terms")
    sqa.termChop(ov2)
    if dbg: 
        for t in ov2:
            print(t)
        print("")

    if dbg: print("remove tags")
    ov_f = []
    for t in ov2:
        tn = sqa.rmTags( t )
        print(tn) 
        ov_f.append( tn )
    if dbg: print("")
    print("")

    return ov_f 


if __name__ == '__main__':
    dbg = True 
    S  = 0
    Ms = 0
    w00_ket  = [] # |MRSF(0,0)> : S=0, Ms=0 MRSF ket vector 
    w00_ket += genCSF([2,2,0,0], S, Ms, ket=True)   # G 
    w00_ket += genCSF([2,1,1,0], S, Ms, ket=True)   # LR 
    w00_ket += genCSF([1,2,1,0], S, Ms, ket=True)   # CO1 
    w00_ket += genCSF([2,1,0,1], S, Ms, ket=True)   # O2V
    w00_ket += genCSF([2,0,2,0], S, Ms, ket=True)   # D 
    w00_ket += genCSF([1,1,2,0], S, Ms, ket=True)   # CO2
    w00_ket += genCSF([2,0,1,1], S, Ms, ket=True)   # O1V
    w00_ket += genCSF([1,1,1,1], S, Ms, ket=True)   # CV

    w00_bra  = [] # <MRSF(0,0)| : S=0, Ms=0 MRSF bra vector 
    w00_bra += genCSF([2,2,0,0], S, Ms, bra=True)   # G O1
    w00_bra += genCSF([2,1,1,0], S, Ms, bra=True)   # LRO1 / O2-O2
    w00_bra += genCSF([1,2,1,0], S, Ms, bra=True)   # CO1 
    w00_bra += genCSF([2,1,0,1], S, Ms, bra=True)   # O2V 
    w00_bra += genCSF([2,0,2,0], S, Ms, bra=True)   # D O2
    w00_bra += genCSF([1,1,2,0], S, Ms, bra=True)   # CO2
    w00_bra += genCSF([2,0,1,1], S, Ms, bra=True)   # O1V
    w00_bra += genCSF([1,1,1,1], S, Ms, bra=True)   # CV

    S  = 1 
    Ms = 0
    w10_ket  = [] # |MRSF(1,0)> : S=1, Ms=0 MRSF ket vector 
    w10_ket += genCSF([2,1,1,0], S, Ms, ket=True)   # O1-O1 / O2-O2
    w10_ket += genCSF([1,2,1,0], S, Ms, ket=True)   # C-O1 
    w10_ket += genCSF([2,1,0,1], S, Ms, ket=True)   # O2-V
    w10_ket += genCSF([1,1,2,0], S, Ms, ket=True)   # C-O2
    w10_ket += genCSF([2,0,1,1], S, Ms, ket=True)   # O1-V
    w10_ket += genCSF([1,1,1,1], S, Ms, ket=True)   # C-V

    w10_bra  = [] # <MRSF(1,0)| : S=1, Ms=0 MRSF bra vector 
    w10_bra += genCSF([2,1,1,0], S, Ms, bra=True)   # O1-O1 / O2-O2
    w10_bra += genCSF([1,2,1,0], S, Ms, bra=True)   # C-O1 
    w10_bra += genCSF([2,1,0,1], S, Ms, bra=True)   # O2-V 
    w10_bra += genCSF([1,1,2,0], S, Ms, bra=True)   # C-O2
    w10_bra += genCSF([2,0,1,1], S, Ms, bra=True)   # O1-V
    w10_bra += genCSF([1,1,1,1], S, Ms, bra=True)   # C-V

    if True:
        print("w00_ket")
        for t in w00_ket:
            print(t)
        print("")
        print("w00_bra")
        for t in w00_bra:
            print(t)
        print("")
        print("w10_ket")
        for t in w10_ket:
            print(t)
        print("")
        print("w10_bra")
        for t in w10_bra:
            print(t)
        print("")

    # overlap of <MRSF(0,0)|MRSF(0,0)> for sanity check
    print("overlap <00|00>")
    ops = []
    S00 = overlap(w00_bra, ops, w00_ket, dbg=dbg)

#    # overlap of <MRSF(1,0)|MRSF(1,0)> for sanity check
#    print "overlap <10|10>"
#    ops = []
#    S11 = overlap(w10_bra, ops, w10_ket, dbg=dbg)
#
#    # overlap of <MRSF(0,0)|MRSF(1,0)> for sanity check
#    print "overlap <00|10>"
#    ops = []
#    S01 = overlap(w00_bra, ops, w10_ket, dbg=dbg)
#
#    # overlap of <MRSF(0,0)|MRSF(1,0)> for sanity check
#    print "overlap <10|00>"
#    ops = []
#    S10 = overlap(w10_bra, ops, w00_ket, dbg=dbg)

#    # spin-dependent transition density matrix 1 
#    # <MRSF(0,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)>
#    print "SDTDM <00|a^+_ua a_ta|10>"
#    ops = [t_ua_cre, t_ta_des]
#    T01aa = overlap(w00_bra, ops, w10_ket, dbg=dbg)
  
#    # spin-dependent transition density matrix 2 
#    # <MRSF(0,0)| a^+_{u beta} a_{t beta} |MRSF(1,0)>
#    print "SDTDM <00|a^+_ub a_tb|10>"
#    ops = [t_ub_cre, t_tb_des]
#    T01bb = overlap(w00_bra, ops, w10_ket, dbg=dbg)

#    # spin-dependent transition density matrix 3 
#    # <MRSF(0,0)| a^+_{u alpha} a_{t alpha} |MRSF(0,0)>
#    print "SDTDM <00|a^+_ua a_ta|00>"
#    ops = [t_ua_cre, t_ta_des]
#    T00aa = overlap(w00_bra, ops, w00_ket, dbg=dbg)
#
#    # spin-dependent transition density matrix 4 
#    # <MRSF(0,0)| a^+_{u beta} a_{t beta} |MRSF(0,0)>
#    print "SDTDM <00|a^+_ub a_tb|00>"
#    ops = [t_ub_cre, t_tb_des]
#    T00bb = overlap(w00_bra, ops, w00_ket, dbg=dbg)
#
#    # spin-dependent transition density matrix 5 
#    # <MRSF(1,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)>
#    print "SDTDM <10|a^+_ua a_ta|10>"
#    ops = [t_ua_cre, t_ta_des]
#    T11aa = overlap(w10_bra, ops, w10_ket, dbg=dbg)
#
#    # spin-dependent transition density matrix 6 
#    # <MRSF(1,0)| a^+_{u beta} a_{t beta} |MRSF(1,0)>
#    print "SDTDM <00|a^+_ub a_tb|00>"
#    ops = [t_ub_cre, t_tb_des]
#    T11bb = overlap(w10_bra, ops, w10_ket, dbg=dbg)
#
