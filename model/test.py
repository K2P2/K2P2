from model import *
model=Model('new')
model.add_susy(0)
su5=LieGroup('su',5)
su3=LieGroup('su',3)
a1=LieGroup('u',1,Rational(1,3))
a2=LieGroup('u',1,Rational(1,2))
ggrp=ProdGroup([su5,su3],[a1,a2])
model.add_ggrp(ggrp)
model.Break(0)
model.Break()
print model.broken_glag[0]

