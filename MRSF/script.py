import secondQuantizationAlgebra as sqa
import math 
import time

sqa.options.verbose = False
# tags 
tag_alpha   = sqa.options.alpha_type
tag_beta    = sqa.options.beta_type

tag_core    = sqa.options.core_type
tag_active1 = sqa.options.active_type
tag_active2 = sqa.options.active_type2
tag_virtual = sqa.options.virtual_type
tag_total   = tag_core + tag_active1 + tag_active2 + tag_virtual

# indices for wavef 
dummy = True
Cca = sqa.index('Cca', [tag_core, tag_alpha],   dummy)
Ccb = sqa.index('Ccb', [tag_core, tag_beta ],   dummy)
Cia = sqa.index('Cia', [tag_core, tag_alpha],   dummy)
Cib = sqa.index('Cib', [tag_core, tag_beta ],   dummy)
Cja = sqa.index('Cja', [tag_core, tag_alpha],   dummy)
Cjb = sqa.index('Cjb', [tag_core, tag_beta ],   dummy)
O1a = sqa.index('O1a', [tag_active1, tag_alpha],dummy)
O1b = sqa.index('O1b', [tag_active1, tag_beta ],dummy)
O2a = sqa.index('O2a', [tag_active2, tag_alpha],dummy)
O2b = sqa.index('O2b', [tag_active2, tag_beta ],dummy)
Vaa = sqa.index('Vaa', [tag_virtual, tag_alpha],dummy)
Vab = sqa.index('Vab', [tag_virtual, tag_beta ],dummy)
Vba = sqa.index('Vba', [tag_virtual, tag_alpha],dummy)
Vbb = sqa.index('Vbb', [tag_virtual, tag_beta ],dummy)

# indices for Op 
dummy = False 
t_a = sqa.index('t_a', [tag_total, tag_alpha],  dummy)
t_b = sqa.index('t_b', [tag_total, tag_beta ],  dummy)
u_a = sqa.index('u_a', [tag_total, tag_alpha],  dummy)
u_b = sqa.index('u_b', [tag_total, tag_beta ],  dummy)

# 
Cca_cre = sqa.creOp(Cca)
Ccb_cre = sqa.creOp(Ccb)
Cja_cre = sqa.creOp(Cja)
Cjb_cre = sqa.creOp(Cjb)
O1a_cre = sqa.creOp(O1a)
O1b_cre = sqa.creOp(O1b)
O2a_cre = sqa.creOp(O2a)
O2b_cre = sqa.creOp(O2b)
Vba_cre = sqa.creOp(Vba)
Vbb_cre = sqa.creOp(Vbb)

Cca_des = sqa.desOp(Cca)
Ccb_des = sqa.desOp(Ccb)
Cia_des = sqa.desOp(Cia)
Cib_des = sqa.desOp(Cib)
O1a_des = sqa.desOp(O1a)
O1b_des = sqa.desOp(O1b)
O2a_des = sqa.desOp(O2a)
O2b_des = sqa.desOp(O2b)
Vaa_des = sqa.desOp(Vaa)
Vab_des = sqa.desOp(Vab)

ta_des = sqa.desOp(t_a)
tb_des = sqa.desOp(t_b)
ua_cre = sqa.creOp(u_a)
ub_cre = sqa.creOp(u_b)

t_ta_des = sqa.term(1.0,[],[ta_des])
t_tb_des = sqa.term(1.0,[],[tb_des])
t_ua_cre = sqa.term(1.0,[],[ua_cre])
t_ub_cre = sqa.term(1.0,[],[ub_cre])


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

  if (S == 0 and Ms == 0) :                                 # Singlet
    #G
    if nc1 == 2 and nc2 == 2 and nv1 == 0 and nv2 == 0:
        CGcoeff = 1.0
        if ket:
            X = [sqa.tensor('X_G         ', [O2a,O1b], [])]
            Op= [O1b_cre, O1a_cre, Ccb_cre, Cca_cre]
            terms = [sqa.term(CGcoeff,[],X+Op)]
        elif bra:
            X = [sqa.tensor('X*_G        ', [O2a,O1b], [])]
            Op= [Cca_des, Ccb_des, O1a_des, O1b_des]
            terms = [sqa.term(CGcoeff,[],X+Op)]
    #D
    elif nc1 == 2 and nc2 == 0 and nv1 == 2 and nv2 == 0:
        CGcoeff = 1.0
        if ket:
            X = [sqa.tensor('X_D         ', [O1b,O2a], [])]
            Op= [O2a_cre, O2b_cre, Ccb_cre, Cca_cre] 
            terms = [sqa.term(CGcoeff,[],X+Op)]
        elif bra:
            X = [sqa.tensor('X*_D        ', [O1b,O2a], [])]
            Op= [Cca_des, Ccb_des, O2b_des, O2a_des]
            terms = [sqa.term(CGcoeff,[],X+Op)]
    #L,R
    elif nc1 == 2 and nc2 == 1 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_LR1       ', [O1b,O1a], [])]
            Op= [O2a_cre, O1b_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_LR2       ', [O2b,O2a], [])]
            Op= [O2b_cre, O1a_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_LR1      ', [O1b,O1a], [])]
            Op= [Cca_des, Ccb_des, O1b_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_LR2      ', [O2b,O2a], [])]
            Op= [Cca_des, Ccb_des, O1a_des, O2b_des]
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))
    #CO1
    elif nc1 == 1 and nc2 == 2 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_jO1 sum_j ', [Cia,O1b], [])]
            Op= [O2a_cre, O1a_cre, Cjb_cre, O1b_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_jO1 sum_j ', [Cia,O1b], [])]
            Op= [O2b_cre, O1b_cre, O1a_cre, Cja_cre]
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_iO1 sum_i', [Cia,O1b], [])]
            Op= [O1b_des, Cib_des, O1a_des, O2a_des] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_iO1 sum_i', [Cia,O1b], [])]
            Op= [Cia_des, O1a_des, O1b_des, O2b_des]
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))
    #CO2
    elif nc1 == 1 and nc2 == 1 and nv1 == 2 and nv2 == 0:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_jO2 sum_j ', [Cia,O2b], [])]
            Op= [O2a_cre, O1a_cre, Cjb_cre, O2b_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_jO2 sum_j ', [Cia,O2b], [])]
            Op= [O2b_cre, O1b_cre, O2a_cre, Cja_cre] 
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_iO2 sum_i', [Cia,O2b], [])]
            Op= [O2b_des, Cib_des, O1a_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_iO2 sum_i', [Cia,O2b], [])]
            Op= [Cia_des, O2a_des, O1b_des, O2b_des]
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))
    #O1V
    elif nc1 == 2 and nc2 == 0 and nv1 == 1 and nv2 == 1:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_O1b sum_b ', [O1a,Vab], [])]
            Op= [O2a_cre, Vbb_cre, Ccb_cre, Cca_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_O1b sum_b ', [O1a,Vab], [])]
            Op= [O2b_cre, Vba_cre, Ccb_cre, Cca_cre] 
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O1a sum_a', [O1a,Vab], [])]
            Op= [Cca_des, Ccb_des, Vab_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_O1a sum_a', [O1a,Vab], [])]
            Op= [Cca_des, Ccb_des, Vaa_des, O2b_des]
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))
    #O2V
    elif nc1 == 2 and nc2 == 1 and nv1 == 0 and nv2 == 1:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_O2b sum_b ', [O2a,Vab], [])]
            Op= [Vbb_cre, O1a_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_O2b sum_b ', [O2a,Vab], [])]
            Op= [Vba_cre, O1b_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O2a sum_a', [O2a,Vab], [])]
            Op= [Cca_des, Ccb_des, O1a_des, Vab_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_O2a sum_a', [O2a,Vab], [])]
            Op= [Cca_des, Ccb_des, O1b_des, Vaa_des]
            terms.append(sqa.term((-1.0)*CGcoeff,[],X+Op))


  elif  (S == 1 and Ms == 0) :                               # Triplet Ms = 0
    #L,R
    if nc1 == 2 and nc2 == 1 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_LR1       ', [O1b,O1a], [])]
            Op= [O2a_cre, O1b_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_LR2       ', [O2b,O2a], [])]
            Op= [O2b_cre, O1a_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_LR1      ', [O1b,O1a], [])]
            Op= [Cca_des, Ccb_des, O1b_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_LR2      ', [O2b,O2a], [])]
            Op= [Cca_des, Ccb_des, O1a_des, O2b_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #CO1
    elif nc1 == 1 and nc2 == 2 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_jO1 sum_j ', [Cia,O1b], [])]
            Op= [O2a_cre, O1a_cre, Cjb_cre, O1b_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_jO1 sum_j ', [Cia,O1b], [])]
            Op= [O2b_cre, O1b_cre, O1a_cre, Cja_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_iO1 sum_i', [Cia,O1b], [])]
            Op= [O1b_des, Cib_des, O1a_des, O2a_des] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_iO1 sum_i', [Cia,O1b], [])]
            Op= [Cia_des, O1a_des, O1b_des, O2b_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #CO2
    elif nc1 == 1 and nc2 == 1 and nv1 == 2 and nv2 == 0:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_jO2 sum_j ', [Cia,O2b], [])]
            Op= [O2a_cre, O1a_cre, Cjb_cre, O2b_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_jO2 sum_j ', [Cia,O2b], [])]
            Op= [O2b_cre, O1b_cre, O2a_cre, Cja_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_iO2 sum_i', [Cia,O2b], [])]
            Op= [O2b_des, Cib_des, O1a_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_iO2 sum_i', [Cia,O2b], [])]
            Op= [Cia_des, O2a_des, O1b_des, O2b_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #O1V
    elif nc1 == 2 and nc2 == 0 and nv1 == 1 and nv2 == 1:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_O1b sum_b ', [O1a,Vab], [])]
            Op= [O2a_cre, Vbb_cre, Ccb_cre, Cca_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_O1b sum_b ', [O1a,Vab], [])]
            Op= [O2b_cre, Vba_cre, Ccb_cre, Cca_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O1a sum_a', [O1a,Vab], [])]
            Op= [Cca_des, Ccb_des, Vab_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_O1a sum_a', [O1a,Vab], [])]
            Op= [Cca_des, Ccb_des, Vaa_des, O2b_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #O2V
    elif nc1 == 2 and nc2 == 1 and nv1 == 0 and nv2 == 1:
        CGcoeff = 1.0/math.sqrt(2.0) 
        terms = []
        if ket:
            X = [sqa.tensor('X_O2b sum_b ', [O2a,Vab], [])]
            Op= [Vbb_cre, O1a_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X_O2b sum_b ', [O2a,Vab], [])]
            Op= [Vba_cre, O1b_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O2a sum_a', [O2a,Vab], [])]
            Op= [Cca_des, Ccb_des, O1a_des, Vab_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
            X = [sqa.tensor('X*_O2a sum_a', [O2a,Vab], [])]
            Op= [Cca_des, Ccb_des, O1b_des, Vaa_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))


  elif (S == 1 and Ms ==+1):                               # Triplet Ms = +1
    #L,R
    if nc1 == 2 and nc2 == 1 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0
        terms = []
        if ket:
            X = [sqa.tensor('X_LR1       ', [O1b,O1a], [])]
            Op= [O2a_cre, O1a_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_LR1      ', [O1b,O1a], [])]
            Op= [Cca_des, Ccb_des, O1a_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #CO1
    elif nc1 == 1 and nc2 == 2 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0
        terms = []
        if ket:
            X = [sqa.tensor('X_jO1 sum_j ', [Cia,O1b], [])]
            Op= [O2a_cre, O1b_cre, O1a_cre, Cja_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_iO1 sum_i', [Cia,O1b], [])]
            Op= [Cia_des, O1a_des, O1b_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #CO2
    elif nc1 == 1 and nc2 == 1 and nv1 == 2 and nv2 == 0:
        CGcoeff = 1.0
        terms = []
        if ket:
            X = [sqa.tensor('X_jO2 sum_j ', [Cia,O2b], [])]
            Op= [O2a_cre, O1a_cre, Cja_cre, O2b_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_iO2 sum_i', [Cia,O2b], [])]
            Op= [O2b_des, Cia_des, O1a_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #O1V
    elif nc1 == 2 and nc2 == 0 and nv1 == 1 and nv2 == 1:
        CGcoeff = 1.0
        terms = []
        if ket:
            X = [sqa.tensor('X_O1b sum_b ', [O1a,Vab], [])]
            Op= [O2a_cre, Vba_cre, Ccb_cre, Cca_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O1a sum_a', [O1a,Vab], [])]
            Op= [Cca_des, Ccb_des, Vaa_des, O2a_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #O2V
    elif nc1 == 2 and nc2 == 1 and nv1 == 0 and nv2 == 1:
        CGcoeff = 1.0
        terms = []
        if ket:
            X = [sqa.tensor('X_O2b sum_b ', [O2a,Vab], [])]
            Op= [Vba_cre, O1a_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O2a sum_a', [O2a,Vab], [])]
            Op= [Cca_des, Ccb_des, O1a_des, Vaa_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))

  elif (S == 1 and Ms ==-1):                               # Triplet Ms = -1
    #L,R
    if nc1 == 2 and nc2 == 1 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0 
        terms = []
        if ket:
            X = [sqa.tensor('X_LR1       ', [O1b,O1a], [])]
            Op= [O2b_cre, O1b_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_LR1      ', [O1b,O1a], [])]
            Op= [Cca_des, Ccb_des, O1b_des, O2b_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #CO1
    elif nc1 == 1 and nc2 == 2 and nv1 == 1 and nv2 == 0:
        CGcoeff = 1.0
        terms = []
        if ket:
            X = [sqa.tensor('X_jO1 sum_j ', [Cia,O1b], [])]
            Op= [O2b_cre, O1a_cre, Cjb_cre, O1b_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_iO1 sum_i', [Cia,O1b], [])]
            Op= [O1b_des, Cib_des, O1a_des, O2b_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #CO2
    elif nc1 == 1 and nc2 == 1 and nv1 == 2 and nv2 == 0:
        CGcoeff = 1.0
        terms = []
        if ket:
            X = [sqa.tensor('X_jO2 sum_j ', [Cia,O2b], [])]
            Op= [O2b_cre, O1b_cre, O2a_cre, Cjb_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_iO2 sum_i', [Cia,O2b], [])]
            Op= [Cib_des, O2a_des, O1b_des, O2b_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #O1V
    elif nc1 == 2 and nc2 == 0 and nv1 == 1 and nv2 == 1:
        CGcoeff = 1.0
        terms = []
        if ket:
            X = [sqa.tensor('X_O1b sum_b ', [O1a,Vab], [])]
            Op= [O2b_cre, Vbb_cre, Ccb_cre, Cca_cre] 
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O1a sum_a', [O1a,Vab], [])]
            Op= [Cca_des, Ccb_des, Vab_des, O2b_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
    #O2V
    elif nc1 == 2 and nc2 == 1 and nv1 == 0 and nv2 == 1:
        CGcoeff = 1.0
        terms = []
        if ket:
            X = [sqa.tensor('X_O2b sum_b ', [O2a,Vab], [])]
            Op= [Vbb_cre, O1b_cre, Ccb_cre, Cca_cre]
            terms.append(sqa.term(CGcoeff,[],X+Op))
        elif bra:
            X = [sqa.tensor('X*_O2a sum_a', [O2a,Vab], [])]
            Op= [Cca_des, Ccb_des, O1b_des, Vab_des]
            terms.append(sqa.term(CGcoeff,[],X+Op))
  return terms 

def overlap(bra, ops, ket, dbg=False):
    startTime = time.time()
    print ""
    print "bra:"
    for t in bra:
        print t
    print "ops:"
    for t in ops:
        print t
    print "ket:"
    for t in ket:
        print t
#    print "<bra|ops|ket>:"

    if dbg: print "expand terms"
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
    if dbg: print " (%.3f seconds)\n" %(time.time()-startTime)

    #    normal order form
    print ""
    if True: print "Making normal order for maximum contraction form, removing delta_{tu}, ..."
    startTime = time.time()
    ov_tmp2= []
    i = 0
    for t in ov_tmp:
        i += 1
        if dbg: print "%d / %d "%(i, len(ov_tmp))
        term = sqa.normalOrder_maxcontraction( t )
        for tn in term:
            if dbg: print tn
        ov_tmp2 += term
    if dbg: print " (%.3f seconds)\n" %(time.time()-startTime)

    #TODO: remove terms by vacuumFermi in sqaSLee
    if dbg: print "remove delta_{tu} and delta"
    if dbg: startTime = time.time()
    ov_tmp3= []
    for t in ov_tmp2:
        tn = sqa.vacuumFermi( t )
        if dbg: print tn 
        ov_tmp3.append( tn )
    if dbg: print ""
    if dbg: print " (%.3f seconds)\n" %(time.time()-startTime)

    if dbg: print "contract form"
    if dbg: startTime = time.time()
    for t in ov_tmp3:
        t.contractDeltaFuncs_mrsf()
        if dbg: print t
    print " (%.3f seconds)\n" %(time.time()-startTime)
    print ""
    
    if dbg: print "remove near zero terms"
    if dbg: startTime = time.time()
    sqa.termChop(ov_tmp3)
    if True: 
        for t in ov_tmp3:
            print t
        print ""
    print ""
    if dbg: print " (%.3f seconds)\n" %(time.time()-startTime)

    if False: 
      print ""
      if True: print "Combine terms"
      sqa.combineTerms(ov_tmp3)
      if True: 
         for t in ov_tmp3:
           print t
      
    return ov_tmp3


if __name__ == '__main__':
    dbg = False
    S  = 0
    Ms = 0
    w00_ket  = [] # |MRSF(0,0)> : S=0, Ms=0 MRSF ket vector 
    w00_ket += genCSF([1,2,1,0], S, Ms, ket=True)   # CO1 
    w00_ket += genCSF([2,1,1,0], S, Ms, ket=True)   # LR1 / LR2
    w00_ket += genCSF([2,2,0,0], S, Ms, ket=True)   # G 
    w00_ket += genCSF([1,1,2,0], S, Ms, ket=True)   # CO2
    w00_ket += genCSF([2,0,2,0], S, Ms, ket=True)   # D 
    w00_ket += genCSF([2,0,1,1], S, Ms, ket=True)   # O1V
    w00_ket += genCSF([2,1,0,1], S, Ms, ket=True)   # O2V
#    w00_ket += genCSF([1,1,1,1], S, Ms, ket=True)   # CV

    w00_bra  = [] # <MRSF(0,0)| : S=0, Ms=0 MRSF bra vector 
    w00_bra += genCSF([1,2,1,0], S, Ms, bra=True)   # CO1 
    w00_bra += genCSF([2,1,1,0], S, Ms, bra=True)   # LR1 / LR2
    w00_bra += genCSF([2,2,0,0], S, Ms, bra=True)   # G O1
    w00_bra += genCSF([1,1,2,0], S, Ms, bra=True)   # CO2
    w00_bra += genCSF([2,0,2,0], S, Ms, bra=True)   # D O2
    w00_bra += genCSF([2,0,1,1], S, Ms, bra=True)   # O1V
    w00_bra += genCSF([2,1,0,1], S, Ms, bra=True)   # O2V 
#    w00_bra += genCSF([1,1,1,1], S, Ms, bra=True)   # CV

    S  = 1 
    Ms = 0
    w10_ket  = [] # |MRSF(1,0)> : S=1, Ms=0 MRSF ket vector 
    w10_ket += genCSF([1,2,1,0], S, Ms, ket=True)   # CO1 
    w10_ket += genCSF([2,1,1,0], S, Ms, ket=True)   # LR1 / LR2
    w10_ket += genCSF([1,1,2,0], S, Ms, ket=True)   # CO2
    w10_ket += genCSF([2,0,1,1], S, Ms, ket=True)   # O1V
    w10_ket += genCSF([2,1,0,1], S, Ms, ket=True)   # O2V
#    w10_ket += genCSF([1,1,1,1], S, Ms, ket=True)   # CV

    w10_bra  = [] # <MRSF(1,0)| : S=1, Ms=0 MRSF bra vector 
    w10_bra += genCSF([1,2,1,0], S, Ms, bra=True)   # CO1 
    w10_bra += genCSF([2,1,1,0], S, Ms, bra=True)   # LR1 / LR2
    w10_bra += genCSF([1,1,2,0], S, Ms, bra=True)   # CO2
    w10_bra += genCSF([2,0,1,1], S, Ms, bra=True)   # O1V
    w10_bra += genCSF([2,1,0,1], S, Ms, bra=True)   # O2V 
#    w10_bra += genCSF([1,1,1,1], S, Ms, bra=True)   # CV

    S  = 1 
    Ms =+1
    w1p1_ket  = [] # |MRSF(1,0)> : S=1, Ms=1 MRSF ket vector 
    w1p1_ket += genCSF([1,2,1,0], S, Ms, ket=True)   # CO1 
    w1p1_ket += genCSF([2,1,1,0], S, Ms, ket=True)   # LR1 / LR2
    w1p1_ket += genCSF([1,1,2,0], S, Ms, ket=True)   # CO2
    w1p1_ket += genCSF([2,0,1,1], S, Ms, ket=True)   # O1V
    w1p1_ket += genCSF([2,1,0,1], S, Ms, ket=True)   # O2V

    w1p1_bra  = [] # <MRSF(1,0)| : S=1, Ms=1 MRSF bra vector 
    w1p1_bra += genCSF([1,2,1,0], S, Ms, bra=True)   # CO1 
    w1p1_bra += genCSF([2,1,1,0], S, Ms, bra=True)   # LR1 / LR2
    w1p1_bra += genCSF([1,1,2,0], S, Ms, bra=True)   # CO2
    w1p1_bra += genCSF([2,0,1,1], S, Ms, bra=True)   # O1V
    w1p1_bra += genCSF([2,1,0,1], S, Ms, bra=True)   # O2V 

    S  = 1 
    Ms =-1
    w1m1_ket  = [] # |MRSF(1,0)> : S=1, Ms=-1 MRSF ket vector 
    w1m1_ket += genCSF([1,2,1,0], S, Ms, ket=True)   # CO1 
    w1m1_ket += genCSF([2,1,1,0], S, Ms, ket=True)   # LR1 / LR2
    w1m1_ket += genCSF([1,1,2,0], S, Ms, ket=True)   # CO2
    w1m1_ket += genCSF([2,0,1,1], S, Ms, ket=True)   # O1V
    w1m1_ket += genCSF([2,1,0,1], S, Ms, ket=True)   # O2V

    w1m1_bra  = [] # <MRSF(1,0)| : S=1, Ms=-1 MRSF bra vector 
    w1m1_bra += genCSF([1,2,1,0], S, Ms, bra=True)   # CO1 
    w1m1_bra += genCSF([2,1,1,0], S, Ms, bra=True)   # LR1 / LR2
    w1m1_bra += genCSF([1,1,2,0], S, Ms, bra=True)   # CO2
    w1m1_bra += genCSF([2,0,1,1], S, Ms, bra=True)   # O1V
    w1m1_bra += genCSF([2,1,0,1], S, Ms, bra=True)   # O2V


    if dbg: 
        print "w00_ket"
        for t in w00_ket:
            print t
        print ""
        print "w00_bra"
        for t in w00_bra:
            print t
        print ""
        print "w10_ket"
        for t in w10_ket:
            print t
        print ""
        print "w10_bra"
        for t in w10_bra:
            print t
        print ""
        print "w1p1_ket"
        for t in w1p1_ket:
            print t
        print ""
        print "w1p1_bra"
        for t in w1p1_bra:
            print t
        print ""

#
#  S00 = <MRSF(0,0) | MRSF(0,0)>
#
#    print "SDTDM <MRSF(0,0) | MRSF(0,0)>"
#    ops = []
#    S00 = overlap(w00_bra, ops, w00_ket, dbg=0)
#
#  T00 = <MRSF(0,0) | MRSF(1,0)> = 0
#
#    print "SDTDM <00|10>"
#    ops = []
#    T00 = overlap(w00_bra, ops, w10_ket, dbg=0)
#
#  T01 = <MRSF(0,0) | MRSF(1,1)> = 0
#
#    print "SDTDM <00|11>"
#    ops = []
#    T01 = overlap(w00_bra, ops, w1p1_ket, dbg=0)
#
#  T0m1 = <MRSF(0,0) | MRSF(1,-1)> = 0
#
#    print "SDTDM <00|1-1>"
#    ops = []
#    T0m1 = overlap(w00_bra, ops, w1m1_ket, dbg=0)
#
#  T110 = <MRSF(1,0) | MRSF(1,0)> 
#
#    print "SDTDM <10|10>"
#    ops = []
#    T110 = overlap(w10_bra, ops, w10_ket, dbg=0)
#
#  T11 = <MRSF(1,0) | MRSF(1,1)> = 0
#
#    print "SDTDM <10|11>"
#    ops = []
#    T11 = overlap(w10_bra, ops, w1p1_ket, dbg=0)
#
#  T1m1 = <MRSF(1,0) | MRSF(1,-1)> = 0
#
#    print "SDTDM <10|1-1>"
#    ops = []
#    T1m1 = overlap(w10_bra, ops, w1m1_ket, dbg=0)
#
#############################################################

#
#  S00aa = <MRSF(0,0)| a^+_{u alpha} a_{t alpha} |MRSF(0,0)>
#
#    print "SDTDM <00|a^+_ua a_ta|00>"
#    ops = [t_ua_cre, t_ta_des]
#    S00aa = overlap(w00_bra, ops, w00_ket, dbg=0)

#
#  S00bb = <MRSF(0,0)| a^+_{u beta} a_{t beta} |MRSF(0,0)>  = S00aa
#
#    print "SDTDM <00|a^+_ub a_tb|00>"
#    ops = [t_ub_cre, t_tb_des]
#    S00bb = overlap(w00_bra, ops, w00_ket, dbg=0)

#
#  S00ab = <MRSF(0,0)| a^+_{u beta} a_{t alpha} |MRSF(0,0)> = 0 
#
#    print "SDTDM <00|a^+_ub a_ta|00>"
#    ops = [t_ub_cre, t_ta_des]
#    S00ab = overlap(w00_bra, ops, w00_ket, dbg=0)

#
#  S00ba = <MRSF(0,0)| a^+_{u alpha} a_{t beta} |MRSF(0,0)> = 0
#
#    print "SDTDM <00|a^+_ua a_tb|00>"
#    ops = [t_ua_cre, t_tb_des]
#    S00ba = overlap(w00_bra, ops, w00_ket, dbg=0)
#############################################################

#
#  T00aa = <MRSF(0,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)> 
#
#    print "SDTDM <00|a^+_ua a_ta|10>"
#    ops = [t_ua_cre, t_ta_des]
#    T00aa = overlap(w00_bra, ops, w10_ket, dbg=0)

#
#  T00bb =  <MRSF(0,0)| a^+_{u beta} a_{t beta} |MRSF(1,0)> = - T00aa
#
#    print "SDTDM <00|a^+_ub a_tb|10>"
#    ops = [t_ub_cre, t_tb_des]
#    T00bb = overlap(w00_bra, ops, w10_ket, dbg=0)

#
#  T00ab = <MRSF(0,0)| a^+_{u beta} a_{t alpha} |MRSF(1,0)> = 0 
#
#    print "SDTDM <00|a^+_ub a_ta|10>"
#    ops = [t_ub_cre, t_ta_des]
#    T00ab = overlap(w00_bra, ops, w10_ket, dbg=0)

#
#  T00ba = <MRSF(0,0)| a^+_{u alpha} a_{t beta} |MRSF(1,0)> = 0
#
#    print "SDTDM <00|a^+_ua a_tb|10>"
#    ops = [t_ua_cre, t_tb_des]
#    T00ba = overlap(w00_bra, ops, w10_ket, dbg=0)
#############################################################

#
#  T01aa = <MRSF(0,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,1)> = 0
#
#    print "SDTDM <00|a^+_ua a_ta|11>"
#    ops = [t_ua_cre, t_ta_des]
#    T01aa = overlap(w00_bra, ops, w1p1_ket, dbg=0)

#
#  T01bb =  <MRSF(0,0)| a^+_{u beta} a_{t beta} |MRSF(1,1)> = 0
#
#    print "SDTDM <00|a^+_ub a_tb|11>"
#    ops = [t_ub_cre, t_tb_des]
#    T01bb = overlap(w00_bra, ops, w1p1_ket, dbg=0)

#
#  T01ab = <MRSF(0,0)| a^+_{u beta} a_{t alpha} |MRSF(1,1)> = -sqrt(2) * T00aa
#
#    print "SDTDM <00|a^+_ub a_ta|11>"
#    ops = [t_ub_cre, t_ta_des]
#    T01ab = overlap(w00_bra, ops, w1p1_ket, dbg=0)

#
#  T01ab_eq_sqrt2_T00aa = - sqrt(2) *  <MRSF(0,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)> 
#
#    t_ua_cre = sqa.term((-1.0)*math.sqrt(2.0),[],[ua_cre])
#    print "SDTDM -sqrt(2) * <00|a^+_ua a_ta|10>"
#    ops = [t_ua_cre, t_ta_des]
#    T01ab_eq_sqrt2_T00aa = overlap(w00_bra, ops, w10_ket, dbg=0)

#
#  T01ba = <MRSF(0,0)| a^+_{u alpha} a_{t beta} |MRSF(1,1)> = 0
#
#    print "SDTDM <00|a^+_ua a_tb|11>"
#    ops = [t_ua_cre, t_tb_des]
#    T01ba = overlap(w00_bra, ops, w1p1_ket, dbg=0)
#############################################################

#
#  T0m1aa = <MRSF(0,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,-1)> = 0 
#
#    print "SDTDM <00|a^+_ua a_ta|1-1>"
#    ops = [t_ua_cre, t_ta_des]
#    T0m1aa = overlap(w00_bra, ops, w1m1_ket, dbg=0)

#
#  T0m1bb =  <MRSF(0,0)| a^+_{u beta} a_{t beta} |MRSF(1,-1)> = 0
#
#    print "SDTDM <00|a^+_ub a_tb|1-1>"
#    ops = [t_ub_cre, t_tb_des]
#    T0m1bb = overlap(w00_bra, ops, w1m1_ket, dbg=0)

#
#  T0m1ab = <MRSF(0,0)| a^+_{u beta} a_{t alpha} |MRSF(1,-1)> = 0
#
#    print "SDTDM <00|a^+_ub a_ta|1-1>"
#    ops = [t_ub_cre, t_ta_des]
#    T0m1ab = overlap(w00_bra, ops, w1m1_ket, dbg=0)

#
#  T0m1ba = <MRSF(0,0)| a^+_{u alpha} a_{t beta} |MRSF(1,-1)> = sqrt(2) * T00aa
#
#    print "SDTDM <00|a^+_ua a_tb|1-1>"
#    ops = [t_ua_cre, t_tb_des]
#    T0m1ba = overlap(w00_bra, ops, w1m1_ket, dbg=0)
#############################################################

#
#  T110aa = <MRSF(1,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)>
#
#    print "SDTDM <10|a^+_ua a_ta|10>"
#    ops = [t_ua_cre, t_ta_des]
#    T110aa = overlap(w10_bra, ops, w10_ket, dbg=0)

#
#  T110bb =  <MRSF(1,0)| a^+_{u beta} a_{t beta} |MRSF(1,0)> = T110aa
#
#    print "SDTDM <10|a^+_ub a_tb|10>"
#    ops = [t_ub_cre, t_tb_des]
#    T110bb = overlap(w10_bra, ops, w10_ket, dbg=0)

#
#  T110ab = <MRSF(1,0)| a^+_{u beta} a_{t alpha} |MRSF(1,0)> = 0
#
#    print "SDTDM <10|a^+_ub a_ta|10>"
#    ops = [t_ub_cre, t_ta_des]
#    T110ab = overlap(w10_bra, ops, w10_ket, dbg=0)

#
#  T110ba = <MRSF(1,0)| a^+_{u alpha} a_{t beta} |MRSF(1,0)> = 0
#
#    print "SDTDM <10|a^+_ua a_tb|10>"
#    ops = [t_ua_cre, t_tb_des]
#    T110ba = overlap(w10_bra, ops, w10_ket, dbg=0)
#############################################################

#
#  T11aa = <MRSF(1,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,1)> = 0
#
#    print "SDTDM <10|a^+_ua a_ta|11>"
#    ops = [t_ua_cre, t_ta_des]
#    T11aa = overlap(w10_bra, ops, w1p1_ket, dbg=0)

#
#  T11bb =  <MRSF(1,0)| a^+_{u beta} a_{t beta} |MRSF(1,1)> = 0
#
#    print "SDTDM <10|a^+_ub a_tb|11>"
#    ops = [t_ub_cre, t_tb_des]
#    T11bb = overlap(w10_bra, ops, w1p1_ket, dbg=0)

#
#  T11ab = <MRSF(1,0)| a^+_{u beta} a_{t alpha} |MRSF(1,1)>
#
#    print "SDTDM <10|a^+_ub a_ta|11>"
#    ops = [t_ub_cre, t_ta_des]
#    T11ab = overlap(w10_bra, ops, w1p1_ket, dbg=0)

#
#  T11ba = <MRSF(1,0)| a^+_{u alpha} a_{t beta} |MRSF(1,1)> = 0
#
#    print "SDTDM <10|a^+_ua a_tb|11>"
#    ops = [t_ua_cre, t_tb_des]
#    T11ba = overlap(w10_bra, ops, w1p1_ket, dbg=0)
#
#############################################################

#
#  T1m1aa = <MRSF(1,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,-1)> = 0
#
#    print "SDTDM <10|a^+_ua a_ta|1-1>"
#    ops = [t_ua_cre, t_ta_des]
#    T1m1aa = overlap(w10_bra, ops, w1m1_ket, dbg=0)

#
#  T1m1bb =  <MRSF(1,0)| a^+_{u beta} a_{t beta} |MRSF(1,-1)> = 0
#
#    print "SDTDM <10|a^+_ub a_tb|1-1>"
#    ops = [t_ub_cre, t_tb_des]
#    T1m1bb = overlap(w10_bra, ops, w1m1_ket, dbg=0)

#
#  T1m1ab = <MRSF(1,0)| a^+_{u beta} a_{t alpha} |MRSF(1,-1)> = 0
#
#    print "SDTDM <10|a^+_ub a_ta|1-1>"
#    ops = [t_ub_cre, t_ta_des]
#    T1m1ab = overlap(w10_bra, ops, w1m1_ket, dbg=0)

#
#  T1m1ba = <MRSF(1,0)| a^+_{u alpha} a_{t beta} |MRSF(1,-1)> = T11ab
#
#    print "SDTDM <10|a^+_ua a_tb|1-1>"
#    ops = [t_ua_cre, t_tb_des]
#    T1m1ba = overlap(w10_bra, ops, w1m1_ket, dbg=0)
#######################################

#
#  T111aa = <MRSF(1,1)| a^+_{u alpha} a_{t alpha} |MRSF(1,1)>  = T110aa + 1/sqrt(2) * T11ab
#
#    print "SDTDM <11|a^+_ua a_ta|11>"
#    ops = [t_ua_cre, t_ta_des]
#    T111aa = overlap(w1p1_bra, ops, w1p1_ket, dbg=0)

#
#  T111aa_eq_T110aa_T11ab = <MRSF(1,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)> + 1/sqrt(2) *  <MRSF(1,0)| a^+_{u beta} a_{t alpha} |MRSF(1,1)> 
#
#    t_ub_cre = sqa.term((1.0)/math.sqrt(2.0),[],[ub_cre])
#    print "SDTDM 1/sqrt(2) * <10|a^+_ub a_ta|11>"
#    ops = [t_ub_cre, t_ta_des]
#    T111aa_eq_T110aa_T11ab = overlap(w10_bra, ops, w1p1_ket, dbg=0)

#
#  T111bb =  <MRSF(1,1)| a^+_{u beta} a_{t beta} |MRSF(1,1)> = T110aa - 1/sqrt(2) * T11ab
#
#    print "SDTDM <11|a^+_ub a_tb|11>"
#    ops = [t_ub_cre, t_tb_des]
#    T111bb = overlap(w1p1_bra, ops, w1p1_ket, dbg=0)

#  T111bb_eq_T110aa_T11ab = <MRSF(1,0)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)> - 1/sqrt(2) *  <MRSF(1,0)| a^+_{u beta} a_{t alpha} |MRSF(1,1)> 
#
#    t_ub_cre = sqa.term((-1.0)/math.sqrt(2.0),[],[ub_cre])
#    print "SDTDM -1/sqrt(2) * <10|a^+_ub a_ta|11>"
#    ops = [t_ub_cre, t_ta_des]
#    T111bb_eq_T110aa_T11ab = overlap(w10_bra, ops, w1p1_ket, dbg=0)
#
#  T111ab = <MRSF(1,1)| a^+_{u beta} a_{t alpha} |MRSF(1,1)> = 0
#
#    print "SDTDM <11|a^+_ub a_ta|11>"
#    ops = [t_ub_cre, t_ta_des]
#    T111ab = overlap(w1p1_bra, ops, w1p1_ket, dbg=0)

#
#  T111ba = <MRSF(1,1)| a^+_{u alpha} a_{t beta} |MRSF(1,1)> = 0
#
#    print "SDTDM <11|a^+_ua a_tb|11>"
#    ops = [t_ua_cre, t_tb_des]
#    T111ba = overlap(w1p1_bra, ops, w1p1_ket, dbg=0)
#############################################################

#
#  T11m1aa = <MRSF(1,1)| a^+_{u alpha} a_{t alpha} |MRSF(1,-1)> = 0
#
#    print "SDTDM <11|a^+_ua a_ta|1-1>"
#    ops = [t_ua_cre, t_ta_des]
#    T11m1aa = overlap(w1p1_bra, ops, w1m1_ket, dbg=0)

#
#  T11m1bb =  <MRSF(1,1)| a^+_{u beta} a_{t beta} |MRSF(1,-1)> = 0
#
#    print "SDTDM <11|a^+_ub a_tb|1-1>"
#    ops = [t_ub_cre, t_tb_des]
#    T11m1bb = overlap(w1p1_bra, ops, w1m1_ket, dbg=0)

#
#  T11m1ab = <MRSF(1,1)| a^+_{u beta} a_{t alpha} |MRSF(1,-1)> = 0
#
#    print "SDTDM <11|a^+_ub a_ta|1-1>"
#    ops = [t_ub_cre, t_ta_des]
#    T11m1ab = overlap(w1p1_bra, ops, w1m1_ket, dbg=0)

#
#  T11m1ba = <MRSF(1,1)| a^+_{u alpha} a_{t beta} |MRSF(1,-1)> = 0
#
#    print "SDTDM <11|a^+_ua a_tb|1-1>"
#    ops = [t_ua_cre, t_tb_des]
#    T11m1ba = overlap(w1p1_bra, ops, w1m1_ket, dbg=0)
#############################################################

#
#  Tm11m1aa = <MRSF(1,-1)| a^+_{u alpha} a_{t alpha} |MRSF(1,-1)> = T110aa - 1/sqrt(2) * T11ab
#
#    print "SDTDM <1-1|a^+_ua a_ta|1-1>"
#    ops = [t_ua_cre, t_ta_des]
#    Tm11m1aa = overlap(w1m1_bra, ops, w1m1_ket, dbg=0)

#
#  Tm11m1bb =  <MRSF(1,-1)| a^+_{u beta} a_{t beta} |MRSF(1,-1)> = T110aa + 1/sqrt(2) * T11ab
#
#    print "SDTDM <1-1|a^+_ub a_tb|1-1>"
#    ops = [t_ub_cre, t_tb_des]
#    Tm11m1bb = overlap(w1m1_bra, ops, w1m1_ket, dbg=0)

#
#  Tm11m1ab = <MRSF(1,-1)| a^+_{u beta} a_{t alpha} |MRSF(1,-1)> = 0
#
#    print "SDTDM <1-1|a^+_ub a_ta|1-1>"
#    ops = [t_ub_cre, t_ta_des]
#    Tm11m1ab = overlap(w1m1_bra, ops, w1m1_ket, dbg=0)

#
#  Tm11m1ba = <MRSF(1,-1)| a^+_{u alpha} a_{t beta} |MRSF(1,-1)> = 0
#
#    print "SDTDM <1-1|a^+_ua a_tb|1-1>"
#    ops = [t_ua_cre, t_tb_des]
#    Tm11m1ba = overlap(w1m1_bra, ops, w1m1_ket, dbg=0)
#############################################################

#
#  T100aa = <MRSF(1,0)| a^+_{u alpha} a_{t alpha} |MRSF(0,0)> = T00aa
#
#    print "SDTDM <10|a^+_ua a_ta|00>"
#    ops = [t_ua_cre, t_ta_des]
#    T100aa = overlap(w10_bra, ops, w00_ket, dbg=0)

#
#  T100bb =  <MRSF(1,0)| a^+_{u beta} a_{t beta} |MRSF(0,0)> = - T00aa
#
#    print "SDTDM <10|a^+_ub a_tb|00>"
#    ops = [t_ub_cre, t_tb_des]
#    T100bb = overlap(w10_bra, ops, w00_ket, dbg=0)

#
#  T100ab = <MRSF(1,0)| a^+_{u beta} a_{t alpha} |MRSF(0,0)> = 0
#
#    print "SDTDM <10|a^+_ub a_ta|00>"
#    ops = [t_ub_cre, t_ta_des]
#    T100ab = overlap(w10_bra, ops, w00_ket, dbg=0)

#
#  T100ba = <MRSF(1,0)| a^+_{u alpha} a_{t beta} |MRSF(0,0)> = 0
#
#    print "SDTDM <10|a^+_ua a_tb|00>"
#    ops = [t_ua_cre, t_tb_des]
#    T100ba = overlap(w10_bra, ops, w00_ket, dbg=0)
#############################################################

#
#  T1100aa = <MRSF(1,1)| a^+_{u alpha} a_{t alpha} |MRSF(0,0)> = 0
#
#    print "SDTDM <11|a^+_ua a_ta|00>"
#    ops = [t_ua_cre, t_ta_des]
#    T1110aa = overlap(w1p1_bra, ops, w00_ket, dbg=0)

#
#  T1100bb = <MRSF(1,1)| a^+_{u beta} a_{t beta} |MRSF(0,0)> = 0
#
#    print "SDTDM <11|a^+_ub a_tb|00>"
#    ops = [t_ub_cre, t_tb_des]
#    T1110bb = overlap(w1p1_bra, ops, w00_ket, dbg=0)

#
#  T1100ab = <MRSF(1,1)| a^+_{u beta} a_{t alpha} |MRSF(0,0)> = 0 
#
#    print "SDTDM <11|a^+_ub a_ta|00>"
#    ops = [t_ub_cre, t_ta_des]
#    T1100ab = overlap(w1p1_bra, ops, w00_ket, dbg=0)

#
#  T1110ba = <MRSF(1,1)| a^+_{u alpha} a_{t beta} |MRSF(0,0)> = -sqrt(2) * T00aa
#
#    print "SDTDM <11|a^+_ua a_tb|00>"
#    ops = [t_ua_cre, t_tb_des]
#    T1100ba = overlap(w1p1_bra, ops, w00_ket, dbg=0)
#############################################################

#
#  T1110aa = <MRSF(1,1)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)> = 0
#
#    print "SDTDM <11|a^+_ua a_ta|10>"
#    ops = [t_ua_cre, t_ta_des]
#    T1110aa = overlap(w1p1_bra, ops, w10_ket, dbg=0)

#
#  T1110bb = <MRSF(1,1)| a^+_{u beta} a_{t beta} |MRSF(1,0)> = 0
#
#    print "SDTDM <11|a^+_ub a_tb|10>"
#    ops = [t_ub_cre, t_tb_des]
#    T1110bb = overlap(w1p1_bra, ops, w10_ket, dbg=0)

#
#  T1110ab = <MRSF(1,1)| a^+_{u beta} a_{t alpha} |MRSF(1,0)> = 0 
#
#    print "SDTDM <11|a^+_ub a_ta|10>"
#    ops = [t_ub_cre, t_ta_des]
#    T1110ab = overlap(w1p1_bra, ops, w10_ket, dbg=0)

#
#  T1110ba = <MRSF(1,1)| a^+_{u alpha} a_{t beta} |MRSF(1,0)> = T11ab 
#
#    print "SDTDM <11|a^+_ua a_tb|10>"
#    ops = [t_ua_cre, t_tb_des]
#    T1110ba = overlap(w1p1_bra, ops, w10_ket, dbg=0)
#############################################################

#
#  T1m100aa = <MRSF(1,-1)| a^+_{u alpha} a_{t alpha} |MRSF(0,0)> = 0
#
#    print "SDTDM <1-1|a^+_ua a_ta|00>"
#    ops = [t_ua_cre, t_ta_des]
#    T1m110aa = overlap(w1m1_bra, ops, w00_ket, dbg=0)

#
#  T1m100bb = <MRSF(1,-1)| a^+_{u beta} a_{t beta} |MRSF(0,0)> = 0
#
#    print "SDTDM <1-1|a^+_ub a_tb|00>"
#    ops = [t_ub_cre, t_tb_des]
#    T1m110bb = overlap(w1m1_bra, ops, w00_ket, dbg=0)

#
#  T1m100ab = <MRSF(1,-1)| a^+_{u beta} a_{t alpha} |MRSF(0,0)> = sqrt(2) * T00aa
#
#   print "SDTDM <1-1|a^+_ub a_ta|00>"
#   ops = [t_ub_cre, t_ta_des]
#   T1m100ab = overlap(w1m1_bra, ops, w00_ket, dbg=0)

#
#  T1m110ba = <MRSF(1,-1)| a^+_{u alpha} a_{t beta} |MRSF(0,0)> = 0
#
#    print "SDTDM <1-1|a^+_ua a_tb|00>"
#    ops = [t_ua_cre, t_tb_des]
#    T1m100ba = overlap(w1m1_bra, ops, w00_ket, dbg=0)
#############################################################

#
#  T1m110aa = <MRSF(1,-1)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)> = 0
#
#    print "SDTDM <1-1|a^+_ua a_ta|10>"
#    ops = [t_ua_cre, t_ta_des]
#    T1m110aa = overlap(w1m1_bra, ops, w10_ket, dbg=0)

#
#  T1m110bb = <MRSF(1,-1)| a^+_{u beta} a_{t beta} |MRSF(1,0)> = 0
#
#    print "SDTDM <1-1|a^+_ub a_tb|10>"
#    ops = [t_ub_cre, t_tb_des]
#    T1m110bb = overlap(w1m1_bra, ops, w10_ket, dbg=0)

#
#  T1m110ab = <MRSF(1,-1)| a^+_{u beta} a_{t alpha} |MRSF(1,0)> = T11ab
#
#    print "SDTDM <1-1|a^+_ub a_ta|10>"
#    ops = [t_ub_cre, t_ta_des]
#    T1m110ab = overlap(w1m1_bra, ops, w10_ket, dbg=0)

#
#  T1m110ba = <MRSF(1,-1)| a^+_{u alpha} a_{t beta} |MRSF(1,0)> = 0
#
#    print "SDTDM <1-1|a^+_ua a_tb|10>"
#    ops = [t_ua_cre, t_tb_des]
#    T1m110ba = overlap(w1m1_bra, ops, w10_ket, dbg=0)
#############################################################

#
#  Tm11aa = <MRSF(1,-1)| a^+_{u alpha} a_{t alpha} |MRSF(1,0)> = 0
#
#    print "SDTDM <1-1|a^+_ua a_ta|10>"
#    ops = [t_ua_cre, t_ta_des]
#    Tm11aa = overlap(w1m1_bra, ops, w10_ket, dbg=0)

#
#  Tm11bb =  <MRSF(1,-1)| a^+_{u beta} a_{t beta} |MRSF(1,0)> = 0
#
#    print "SDTDM <1-1|a^+_ub a_tb|10>"
#    ops = [t_ub_cre, t_tb_des]
#    Tm11bb = overlap(w1m1_bra, ops, w10_ket, dbg=0)

#
#  Tm11ab = <MRSF(1,-1)| a^+_{u beta} a_{t alpha} |MRSF(1,0)> = T11ab
#
#    print "SDTDM <1-1|a^+_ub a_ta|10>"
#    ops = [t_ub_cre, t_ta_des]
#    Tm11ab = overlap(w1m1_bra, ops, w10_ket, dbg=0)

#
#  Tm11ba = <MRSF(1,-1)| a^+_{u alpha} a_{t beta} |MRSF(1,0)> = 0
#
#    print "SDTDM <1-1|a^+_ua a_tb|10>"
#    ops = [t_ua_cre, t_tb_des]
#    Tm11ba = overlap(w1m1_bra, ops, w10_ket, dbg=0)
#############################################################
