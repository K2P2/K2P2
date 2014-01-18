from model import *
model=Model('new')
model.add_susy(0)
su3=LieGroup('su',3)
g2=LieGroup('g',2)
a1=LieGroup('u',1,Rational(1,3))
a2=LieGroup('u',1,Rational(1,2))
ggrp=ProdGroup([su3,g2],[a1,a2])
model.add_ggrp(ggrp)

