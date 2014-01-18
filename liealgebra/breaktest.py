from user import *
grp=LieGroup('su',6)
rep=Rep(grp,grp.find_adjoint())
mat=rep.mat()
for key in mat:
	print key, (mat[key].H*mat[key]).trace()

