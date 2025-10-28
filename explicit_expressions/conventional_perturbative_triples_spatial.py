import secondQuantizationAlgebra as sqa

sqa.options.verbose = False
# definitions
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_t = tg_c + tg_a + tg_v

i  = [sqa.index('i%i' %j, [tg_c], True) for j in range(8)]
r  = [sqa.index('r%i' %j, [tg_a], True) for j in range(8)]
a  = [sqa.index('a%i' %j, [tg_v], True) for j in range(8)]
m  = [sqa.index('m%i' %j, [tg_t], True) for j in range(8)]

t_sym = [ sqa.symmetry((1,0,3,2),1) ]
W_sym = [ sqa.symmetry((1,0,3,2),1), sqa.symmetry((0,3,2,1),1) ]

tW_ts_ket_1 = [ sqa.tensor('t2Wk', [a[0], a[1], a[2], i[0], i[1], i[2]], []) ]
tW_ops_ket_1= [ sqa.sfExOp([i[0], i[1], i[2], a[0], a[1], a[2]]) ] 

tW_ts_bra_1 = [ sqa.tensor('t2Wb', [i[3], i[4], i[5], a[3], a[4], a[5]], []) ]
tW_ops_bra_1= [ sqa.sfExOp([a[3], a[4], a[5], i[3], i[4], i[5]]) ] 

#
#t_ts_ket_2 = [ sqa.tensor('t2', [a[0], a[1], i[0], i[3]], []) ]
#t_ops_ket_2= [ sqa.sfExOp([i[0], i[3], a[0], a[1]]) ] 
#W_ts_ket_2 = [ sqa.tensor('W', [i[3], a[2], i[1], i[2]], W_sym) ]
#W_ops_ket_2= [ sqa.sfExOp([i[1], i[2], i[3], a[2]]) ] 
#
#t_ts_bra_1 = [ sqa.tensor('t2', [i[4], i[5], a[4], a[7]], t_sym) ]
#t_ops_bra_1= [ sqa.sfExOp([a[4], a[7], i[4], i[5]]) ] 
#W_ts_bra_1 = [ sqa.tensor('W', [a[7], i[6], a[5], a[6]], W_sym) ]
#W_ops_bra_1= [ sqa.sfExOp([a[5], a[6], a[7], i[6]]) ] 
#
#t_ts_bra_2 = [ sqa.tensor('t2', [i[4], i[7], a[4], a[5]], t_sym) ]
#t_ops_bra_2= [ sqa.sfExOp([a[4], a[5], i[4], i[7]]) ] 
#W_ts_bra_2 = [ sqa.tensor('W', [i[5], i[6], i[7], a[6]], W_sym) ]
#W_ops_bra_2= [ sqa.sfExOp([i[7], a[6], i[5], i[6]]) ] 

terms = []
Op1   = tW_ts_ket_1 + tW_ops_ket_1 + tW_ts_bra_1 + tW_ops_bra_1
terms.append( sqa.term( 1.0, [], Op1) )
#Op2 = t_ts_ket_2 + t_ops_ket_2 + W_ts_ket_2 + W_ops_ket_2 + W_ts_bra_1 + W_ops_bra_1 + t_ts_bra_1 + t_ops_bra_1
#terms.append( sqa.term(-1.0, [], Op2) )
#Op3 = t_ts_ket_1 + t_ops_ket_1 + W_ts_ket_1 + W_ops_ket_1 + W_ts_bra_2 + W_ops_bra_2 + t_ts_bra_2 + t_ops_bra_2
#terms.append( sqa.term(-1.0, [], Op3) )
#Op4 = t_ts_ket_2 + t_ops_ket_2 + W_ts_ket_2 + W_ops_ket_2 + W_ts_bra_2 + W_ops_bra_2 + t_ts_bra_2 + t_ops_bra_2
#terms.append( sqa.term( 1.0, [], Op4) )

print("pure")
for t in terms:
    print(t)
print("")

print("normal order form")
Nterms = []
for t in terms:
    Nterm = sqa.normalOrder( t )
    for tn in Nterm:
        print(tn)
    Nterms += Nterm
print("")

print("re-obtain spin-free operator")
Nterms2 = []
for t in Nterms:
    tmp = sqa.SpinFree( t )
    Nterms2.append( tmp )
    print(tmp, tmp.tensors[-1].spin) 
print("")


print("HF Fermi expectation form")
result = []
for t in Nterms2:
    tmp = sqa.HFFermi( t )
    result.append( tmp )
    print(tmp)
print("")

print("contract form")
for t in result:
    t.contractDeltaFuncs()
    print(t)
print("")

print("remove near zero terms")
sqa.termChop(result)
for t in result:
    print(t)
print("")

#sqa.combineTerms(result, maxThreads = 5)
sqa.combineTerms(result)
print("combined terms")
for t in result:
    print(t)
print("")

print("change idx")
result2 = []
for t in result:
    tmp = sqa.idxname( t )
    result2.append( tmp )
    print(tmp)
print("")

#print "obtain weight factors"
#result3 = []
#for t in result2:
#    tmp = sqa.weight_factor( t )
#    result3.append( tmp )
#    print tmp
#print ""

print("non eq idx sort")
result = sqa.sort_noeqidx_terms ( result2 )
for t in result:
    print(t)
print("")


