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
n  = [sqa.index('n%i' %j, [tg_t], True) for j in range(8)]

t_sym = [ sqa.symmetry((1,0,3,2),1) ]
W_sym = [ sqa.symmetry((1,0,3,2),1), sqa.symmetry((0,3,2,1),1) ]

W_ts_ket = [ sqa.tensor('W', [a[2], a[3], i[2], i[3]], W_sym) ]
W_ops_ket= [ sqa.sfExOp([i[2], i[3], a[2], a[3]]) ] 

t_ts_bra = [ sqa.tensor('t2', [i[0], i[1], a[0], a[1]], t_sym) ]
t_ops_bra= [ sqa.sfExOp([a[0], a[1], i[0], i[1]]) ] 

terms = []
#Op1   = t_ts_ket_1 + t_ops_ket_1 + W_ts_ket_1 + W_ops_ket_1 + W_ts_bra_1 + W_ops_bra_1 + t_ts_bra_1 + t_ops_bra_1
#terms.append( sqa.term( 1.0, [], Op1) )
#Op2 = t_ts_ket_2 + t_ops_ket_2 + W_ts_ket_2 + W_ops_ket_2 + W_ts_bra_1 + W_ops_bra_1 + t_ts_bra_1 + t_ops_bra_1
#terms.append( sqa.term( 1.0, [], Op2) )
#Op3 = t_ts_ket_1 + t_ops_ket_1 + W_ts_ket_1 + W_ops_ket_1 + W_ts_bra_2 + W_ops_bra_2 + t_ts_bra_2 + t_ops_bra_2
#terms.append( sqa.term( 1.0, [], Op3) )
#Op4 = t_ts_ket_2 + t_ops_ket_2 + W_ts_ket_2 + W_ops_ket_2 + W_ts_bra_2 + W_ops_bra_2 + t_ts_bra_2 + t_ops_bra_2
#Op4 = t_ops_ket_2 + W_ops_ket_2 + W_ops_bra_2 + t_ops_bra_2
#Op4 = t_ops_ket_2 + W_ops_ket_2 
Op4 = W_ts_ket + W_ops_ket + t_ts_bra + t_ops_bra
#Op4 = W_ops_bra_4 + t_ops_bra_4
terms.append( sqa.term( 1.0, [], Op4) )

print "pure"
for t in terms:
    print t
print ""

print "normal order form"
Nterms = []
for t in terms:
    Nterm = sqa.normalOrder( t )
    for tn in Nterm:
        print tn, tn.tensors[-1].spin 
    Nterms += Nterm
print ""

print "re-obtain spin-free operator"
Nterms2 = []
for t in Nterms:
    tmp = sqa.SpinFree( t )
    Nterms2.append( tmp )
    print tmp, tmp.tensors[-1].spin 
print ""

print "HF Fermi expectation form"
result = []
for t in Nterms2:
    tmp = sqa.HFFermi( t )
    result.append( tmp )
    print tmp
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

#sqa.combineTerms(result, maxThreads = 5)
#sqa.combineTerms(result)
#print "combined terms"
#for t in result:
#    print t
#print ""

print "change idx"
result2 = []
for t in result:
    tmp = sqa.idxname( t )
    result2.append( tmp )
    print tmp
print ""

#print "obtain weight factors"
#result3 = []
#for t in result2:
#    tmp = sqa.weight_factor( t )
#    result3.append( tmp )
#    print tmp
#print ""

print "non eq idx sort"
result = sqa.sort_noeqidx_terms ( result2 )
for t in result:
    print t
print ""


