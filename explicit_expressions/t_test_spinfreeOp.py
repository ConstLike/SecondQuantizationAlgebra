import secondQuantizationAlgebra as sqa

sqa.options.verbose = False
# definitions
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_t = tg_c + tg_a + tg_v

i  = [sqa.index('i%i' %j, [tg_c], True) for j in range(4)]
r  = [sqa.index('r%i' %j, [tg_a], True) for j in range(4)]
a  = [sqa.index('a%i' %j, [tg_v], True) for j in range(4)]
m  = [sqa.index('m%i' %j, [tg_t], True) for j in range(8)]

t_sym = [ sqa.symmetry((1,0,3,2),1) ]
W_sym = [ sqa.symmetry((1,0,3,2),1), sqa.symmetry((0,3,2,1),1) ]

t_ts1 = [ sqa.tensor('t2', [a[0], a[1], i[0], i[1]], t_sym) ]
t_ops1= [ sqa.sfExOp([i[0], i[1], a[0], a[1]]) ] 

W_ts1 = [ sqa.tensor('W', m[0:4], W_sym) ]
W_ops1= [ sqa.sfExOp([m[2], m[3], m[1], m[0]]) ] 

W_ts2 = [ sqa.tensor('W', m[4:8], W_sym) ]
W_ops2= [ sqa.sfExOp([m[6], m[7], m[5], m[4]]) ]

t_ts2 = [ sqa.tensor('t2', [i[2], i[3], a[2], a[3]], t_sym) ]
t_ops2= [ sqa.sfExOp([a[2], a[3], i[2], i[3]]) ] 

Op = t_ts1 + t_ops1 + W_ts1 + W_ops1 + W_ts2 + W_ops2 + t_ts2 + t_ops2

term = sqa.term(1.0, [], Op)
print("pure")
print(term)
print("")

print("normal order form")
Nterm = sqa.normalOrder( term )
for t in Nterm:
    print(t)
print("")

print("HF Fermi expectation form")
result = []
for t in Nterm:
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

