from user import *
grp=LieGroup('so',10)
s=grp.sroots()
print s
rep=Rep(grp,grp.find_adjoint())
w=rep.weights()


m=rep.mat()
 

def col(n):
	a=Matrix(rep.dim(),1,zeroes)
	a[n]=1
	return a
def zeroes(i,j):
	return 0


lm=[]
lmi=[]
for i in range(grp.dim):
	tup=tuple(0 for j in range(i))+(1,)+tuple(0 for j in range(grp.dim-i-1))
	lm.append(m[tup])
	r=tup
	res=tuple(0 for j in range(grp.dim))
	for j in range(len(r)):
		res=add(res,mult(s[j],r[j]))
	ad=[0 for j in range(grp.dim)]
	for j in range(grp.dim):
		ad[j]=2*vector_product(res,s[j])/vector_product(s[j],s[j])
	for j in range(len(w)):
		if w[j]==tuple(ad):
			lmi.append(j)

lwi=[]
for i in range(len(w)):
	if w[i]==tuple(0 for i in range(grp.dim)):
		lwi.append(i)

vec=[]

check=0

def normalize(vect):
	norm=0
	for i in range(len(vect)):
		norm+=vect[i]**2
	return ((1/sqrt(norm))*vect,sqrt(norm))

def find_nz(vect):
	for i in range(len(vect)):
		if vect[i]!=0:
			return vect[i]

simple=[]

for i in range(grp.dim):
	simple.append([0 for k in range(grp.dim)])
	for j in range(grp.dim):
		if (normalize(lm[i]*col(lwi[j]))[0]==col(lmi[i])) or (normalize(lm[i]*col(lwi[j]))[0]==-col(lmi[i])):
			simple[-1][j]=-find_nz(lm[i]*col(lwi[j]))
			print i,j
		elif normalize(lm[i]*col(lwi[j]))[1]==0:
			simple[-1][j]=normalize(lm[i]*col(lwi[j]))[1]
			print i,j
		else:
			print lm[i]*col(lwi[j]),col(lmi[i])
			print 'no'
			check=1

if check==0:
	print 'gut!'
	print simple
	for i in range(len(simple)):
		print 'length of %d is' % i, vector_product(simple[i],simple[i])
		for j in range(len(simple)):
			if i>j:
				print i,j,vector_product(simple[i],simple[j])/sqrt(vector_product(simple[i],simple[i])*vector_product(simple[j],simple[j]))

#print m


