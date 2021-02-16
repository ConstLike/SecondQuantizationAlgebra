import secondQuantizationAlgebra as sqa

sqa.options.verbose = False
# definitions
tg_alpha = sqa.options.alpha_type
tg_beta  = sqa.options.beta_type
tg_spin  = tg_alpha + tg_beta

tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_t = tg_c + tg_a + tg_v

ia = [sqa.index('i%i' %j, [tg_c, tg_alpha], True) for j in range(4)]
ib = [sqa.index('i%i' %j, [tg_c, tg_beta ], True) for j in range(4)]
i  = [sqa.index('i%i' %j, [tg_c, tg_spin ], True) for j in range(4)]
ra = [sqa.index('r%i' %j, [tg_a, tg_alpha], True) for j in range(4)]
rb = [sqa.index('r%i' %j, [tg_a, tg_beta ], True) for j in range(4)]
r  = [sqa.index('r%i' %j, [tg_a, tg_spin ], True) for j in range(4)]
aa = [sqa.index('a%i' %j, [tg_v, tg_alpha], True) for j in range(4)]
ab = [sqa.index('a%i' %j, [tg_v, tg_beta ], True) for j in range(4)]
a  = [sqa.index('a%i' %j, [tg_v, tg_spin ], True) for j in range(4)]
ma = [sqa.index('m%i' %j, [tg_t, tg_alpha], True) for j in range(8)]
mb = [sqa.index('m%i' %j, [tg_t, tg_beta ], True) for j in range(8)]
m  = [sqa.index('m%i' %j, [tg_t, tg_spin ], True) for j in range(8)]

t_sym = [ sqa.symmetry((1,0,3,2),1) ]
W_sym = [ sqa.symmetry((1,0,3,2),1), sqa.symmetry((0,3,2,1),1) ]

taaaa_ts1 = [ sqa.tensor('t2', [aa[0], aa[1], ia[0], ia[1]], t_sym) ]
tabab_ts1 = [ sqa.tensor('t2', [aa[0], ab[1], ia[0], ib[1]], t_sym) ]
tbaba_ts1 = [ sqa.tensor('t2', [ab[0], aa[1], ib[0], ia[1]], t_sym) ]
tbbbb_ts1 = [ sqa.tensor('t2', [ab[0], ab[1], ib[0], ib[1]], t_sym) ]

taaaa_ops1= [ sqa.creOp(ia[0]), sqa.creOp(ia[1]), sqa.desOp(aa[1]), sqa.desOp(aa[0]) ] 
tabab_ops1= [ sqa.creOp(ia[0]), sqa.creOp(ib[1]), sqa.desOp(ab[1]), sqa.desOp(aa[0]) ] 
tbaba_ops1= [ sqa.creOp(ib[0]), sqa.creOp(ia[1]), sqa.desOp(aa[1]), sqa.desOp(ab[0]) ] 
tbbbb_ops1= [ sqa.creOp(ib[0]), sqa.creOp(ib[1]), sqa.desOp(ab[1]), sqa.desOp(ab[0]) ] 

W_ts1 = [ sqa.tensor('W', m[0:4], W_sym) ]
W_ops1= [ sqa.creOp(m[0]), sqa.creOp(m[1]), sqa.desOp(m[3]), sqa.desOp(m[2]) ] 

W_ts2 = [ sqa.tensor('W', m[4:8], W_sym) ]
W_ops2= [ sqa.creOp(m[4]), sqa.creOp(m[5]), sqa.desOp(m[7]), sqa.desOp(m[6]) ] 

taaaa_ts2 = [ sqa.tensor('t2', [ia[2], ia[3], aa[2], aa[3]], t_sym) ]
tabab_ts2 = [ sqa.tensor('t2', [ia[2], ib[3], aa[2], ab[3]], t_sym) ]
tbaba_ts2 = [ sqa.tensor('t2', [ib[2], ia[3], ab[2], aa[3]], t_sym) ]
tbbbb_ts2 = [ sqa.tensor('t2', [ib[2], ib[3], ab[2], ab[3]], t_sym) ]

taaaa_ops2= [ sqa.creOp(aa[2]), sqa.creOp(aa[3]), sqa.desOp(ia[3]), sqa.desOp(ia[2]) ] 
tabab_ops2= [ sqa.creOp(aa[2]), sqa.creOp(ab[3]), sqa.desOp(ib[3]), sqa.desOp(ia[2]) ] 
tbaba_ops2= [ sqa.creOp(ab[2]), sqa.creOp(aa[3]), sqa.desOp(ia[3]), sqa.desOp(ib[2]) ] 
tbbbb_ops2= [ sqa.creOp(ab[2]), sqa.creOp(ab[3]), sqa.desOp(ib[3]), sqa.desOp(ib[2]) ] 

aaaaOp = taaaa_ts1 + taaaa_ops1 + W_ts1 + W_ops1 + W_ts2 + W_ops2 + taaaa_ts2 + taaaa_ops2
ababOp = tabab_ts1 + tabab_ops1 + W_ts1 + W_ops1 + W_ts2 + W_ops2 + tabab_ts2 + tabab_ops2
babaOp = tbaba_ts1 + tbaba_ops1 + W_ts1 + W_ops1 + W_ts2 + W_ops2 + tbaba_ts2 + tbaba_ops2
bbbbOp = tbbbb_ts1 + tbbbb_ops1 + W_ts1 + W_ops1 + W_ts2 + W_ops2 + tbbbb_ts2 + tbbbb_ops2

term_aaaa = sqa.term(1.0, [], aaaaOp)
print "pure"
print term_aaaa
print ""

print "normal order form"
result = sqa.normalOrder( term_aaaa )
for t in result:
    print t
print ""

print "contract form"
for t in result:
    t.contractDeltaFuncs()
    print t
print ""

print "remove near zero terms"
sqa.termChop(result)
for t in result:
    print t
print ""

sqa.combineTerms(result)
print "combined terms"
for t in result:
    print t
print ""


