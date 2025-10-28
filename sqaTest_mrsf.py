#!/usr/bin/env python3
"""
MRSF T11ab Benchmark Test

This test reproduces the T11ab calculation from the MRSF script:
T11ab = <MRSF(1,0)| a^+_{u beta} a_{t alpha} |MRSF(1,1)>

Expected output: 22 terms, runtime ~2.3 seconds
"""

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
    print("")
    print("bra:")
    for t in bra:
        print(t)
    print("ops:")
    for t in ops:
        print(t)
    print("ket:")
    for t in ket:
        print(t)

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
    if dbg: print(" (%.3f seconds)\n" %(time.time()-startTime))

    print("")
    if True: print("Making normal order for maximum contraction form, removing delta_{tu}, ...")
    startTime = time.time()
    ov_tmp2= []
    i = 0
    for t in ov_tmp:
        i += 1
        if dbg: print("%d / %d "%(i, len(ov_tmp)))
        term = sqa.normalOrder_maxcontraction( t )
        for tn in term:
            if dbg: print(tn)
        ov_tmp2 += term
    if dbg: print(" (%.3f seconds)\n" %(time.time()-startTime))

    if dbg: print("remove delta_{tu} and delta")
    if dbg: startTime = time.time()
    ov_tmp3= []
    for t in ov_tmp2:
        tn = sqa.vacuumFermi( t )
        if dbg: print(tn)
        ov_tmp3.append( tn )
    if dbg: print("")
    if dbg: print(" (%.3f seconds)\n" %(time.time()-startTime))

    if dbg: print("contract form")
    if dbg: startTime = time.time()
    for t in ov_tmp3:
        t.contractDeltaFuncs_mrsf()
        if dbg: print(t)
    print(" (%.3f seconds)\n" %(time.time()-startTime))
    print("")

    if dbg: print("remove near zero terms")
    if dbg: startTime = time.time()
    sqa.termChop(ov_tmp3)
    if True:
        for t in ov_tmp3:
            print(t)
        print("")
    print("")
    if dbg: print(" (%.3f seconds)\n" %(time.time()-startTime))

    if False:
      print("")
      if True: print("Combine terms")
      sqa.combineTerms(ov_tmp3)
      if True:
         for t in ov_tmp3:
           print(t)

    return ov_tmp3


if __name__ == '__main__':
    print("="*70)
    print("MRSF T11ab Benchmark Test")
    print("="*70)

    total_start = time.time()

    # Generate MRSF(1,0) bra state
    S  = 1
    Ms = 0
    w10_bra  = []
    w10_bra += genCSF([1,2,1,0], S, Ms, bra=True)   # CO1
    w10_bra += genCSF([2,1,1,0], S, Ms, bra=True)   # LR1 / LR2
    w10_bra += genCSF([1,1,2,0], S, Ms, bra=True)   # CO2
    w10_bra += genCSF([2,0,1,1], S, Ms, bra=True)   # O1V
    w10_bra += genCSF([2,1,0,1], S, Ms, bra=True)   # O2V

    # Generate MRSF(1,+1) ket state
    S  = 1
    Ms = +1
    w1p1_ket  = []
    w1p1_ket += genCSF([1,2,1,0], S, Ms, ket=True)   # CO1
    w1p1_ket += genCSF([2,1,1,0], S, Ms, ket=True)   # LR1 / LR2
    w1p1_ket += genCSF([1,1,2,0], S, Ms, ket=True)   # CO2
    w1p1_ket += genCSF([2,0,1,1], S, Ms, ket=True)   # O1V
    w1p1_ket += genCSF([2,1,0,1], S, Ms, ket=True)   # O2V

    # T11ab = <MRSF(1,0)| a^+_{u beta} a_{t alpha} |MRSF(1,1)>
    print("\nSDTDM <10|a^+_ub a_ta|11>")
    ops = [t_ub_cre, t_ta_des]
    T11ab = overlap(w10_bra, ops, w1p1_ket, dbg=0)

    total_time = time.time() - total_start

    print("="*70)
    print(f"Total execution time: {total_time:.3f} seconds")
    print(f"Number of terms: {len(T11ab)}")
    print("="*70)

    if len(T11ab) == 22:
        print("\nTest PASSED: Found expected 22 terms")
    else:
        print(f"\nTest WARNING: Expected 22 terms, found {len(T11ab)}")
