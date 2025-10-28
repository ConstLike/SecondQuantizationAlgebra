import secondQuantizationAlgebra as sqa

tg_a = sqa.options.alpha_type
tg_b = sqa.options.beta_type
tg_c = sqa.options.core_type
tg_v = sqa.options.virtual_type
tg_h = sqa.options.core_type + sqa.options.virtual_type
p0 = sqa.index('p0', [tg_h], True)
p1 = sqa.index('p1', [tg_h], True)
p2 = sqa.index('p2', [tg_h], True)
p3 = sqa.index('p3', [tg_h], True)
q0 = sqa.index('q0', [tg_v], True)
q1 = sqa.index('q1', [tg_c], True)
sym_h2abab = sqa.symmetry((2,3,0,1),1)
h2abab = sqa.tensor('h2abab', [p0, p1, p2, p3], [sym_h2abab])
p0_c = sqa.creOp(p0)
p1_c = sqa.creOp(p1)
p2_d = sqa.desOp(p2)
p3_d = sqa.desOp(p3)
t1aa = sqa.tensor('t1aa', [q0, q1], [])
q0_c = sqa.creOp(q0)
q1_d = sqa.desOp(q1)
term = sqa.term(1.0, [], [t1aa, q0_c, q1_d, h2abab, p0_c, p1_c, p3_d, p2_d])
print(term)
print('')
result = sqa.normalOrder(term)
for t in result:
  print(t)
print('')
for t in result:
  t.contractDeltaFuncs()
  print(t)
print('')
print(result[1].tensors[1].indices[0].indType)
sqa.termChop(result)
for t in result:
  print(t)
print('')
sqa.combineTerms(result)
print("combined terms")
for t in result:
    print(t)
print("")
print(result[1].tensors[1].indices[0].indType)
print("")
print("change idx")
result2 = []
for t in result:
    tmp = sqa.idxname( t )
    result2.append( tmp )
    print(tmp)
print("")


