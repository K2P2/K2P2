from sympy import *
from networkx import *


#-------------------------secondary functions (ACCESSED internally BY PRIMARY FUNCTIONS)-----------------------# 

lines_dict={-Rational(1,2):1,-1/sqrt(2):2,-sqrt(3)/2:3,Rational(0,1):0}

def step(x):
	if x==0:
		return 0
	else:
		return 1

def subtract(a,b):
	if type(a)==list:
		l=[]
		for i in range(len(a)):
			l.append(a[i]-b[i])
		return l
	elif type(a)==tuple:
		l=()
		for i in range(len(a)):	
			l+=(a[i]-b[i],)
		return l

def add(a,b):
	if type(a)==list:
		l=[]
		for i in range(len(a)):
			l.append(a[i]+b[i])
		return l
	elif type(a)==tuple:
		l=()
		for i in range(len(a)):	
			l+=(a[i]+b[i],)
		return l

def mult(List, number):
	l=[]
	for i in range(len(List)):
		l.append(number*List[i])
	return l

def even(x):
	if int(x/2.)==x/2.:
		return True
	else:
		return False

def increase(a,b):
	return a[0:b-1]+(a[b-1]+1,)+a[b:]

def vector_product(a,b):
	prod=0
	for i in range(len(a)):
		prod+=a[i]*b[i]
	return prod

def check1(a,b,c,n):
	vector1=vector[b-1]
	prod=0
	for i in range(len(a)):
		vector2=vector[i]
		prod+=a[i]*vector_product(vector1,vector2)
	q=c[b-1]
	prod=(-2*prod)/vector_product(vector1,vector1)+q
	if n<=total:
		if abs(round(prod)-prod)>1.e-3:
			raise Exception, "The simple roots are not part of a pi-system as the SU(2) index of simple root no %d with respect to simple root no %d is not an integer" % n,b 
	return prod

def calculate_lines(a,b):
	vector1=vector[b-1]
	lowest_vector=[0 for i in range(total)]
	for i in range(len(a)):
		lowest_vector=add(lowest_vector,mult(vector[i],a[i])
)
	prod=0
	for i in range(len(a)):
		vector2=vector[i]
		prod+=a[i]*vector_product(vector1,vector2)
	prod=prod/sqrt(vector_product(vector1,vector1)*vector_product(lowest_vector,lowest_vector))
	return lines_dict[prod]

def calculate_length(a):
	lowest_vector=[0 for i in range(total)]
	for i in range(len(a)):
		lowest_vector=add(lowest_vector,mult(vector[i],a[i])
)
	return sqrt(vector_product(lowest_vector,lowest_vector))

def recur(l):
	coeff=l[0]
	a=l[1]
	b=l[2]
	if b>total:
		add_list=com_dict[b][1]
		coeff1=com_dict[b][0]
		return [coeff*coeff1,[a,add_list]]
	else:
		return [l[0],[l[1],l[2]]]

def root_mat_func(i,j):
	return 2*vector[i][j]/vector_product(vector[i],vector[i])

def result_mat_func(i,j):
	if i==now:	
		return 1
	else:	
		return 0


def delete_last(var):
	for element in prim:
		prim[element][4]=[]
	for element in mod:
		mod[element][4]=[]
	if var=='high':
		for element in weight_dict_high:
			weight_dict_high[element][4]=[]
	elif var=='low':
		for element in weight_dict_low:
			weight_dict_low[element][4]=[]


def scalar_product(a,b):
	#print scalar_product_dict
	if a==[] and b==[]:
		#print "ba"
		return 1
	elif len(a)!=len(b):
		#print "in"
		return 0
	else:
		try:
			a_tup=tuple(a[i] for i in range(len(a)))
			b_tup=tuple(b[i] for i in range(len(b)))
			#print scalar_product_dict[a_tup][b_tup]
			return scalar_product_dict[a_tup][b_tup]
		except KeyError:
			result=0
			#print a
			a_copy=a[:]
			del a_copy[0]
			#print a_copy,a
			marker=0
			for i in range(len(b)):
				if b[i]==a[0]:
					#print "in"
					new_vector=highest_weight_vector
					for j in range(i+1,len(b)):
						new_vector=subtract(new_vector,vector[b[j]-1])
					product=vector_product(new_vector,vector[a[0]-1])
					b_copy=[b[k] for k in range(i)]+[b[k] for k in range(i+1,len(b))]
					result+=product*scalar_product(a_copy,b_copy)
					marker=1
			#print result
			try:
				scalar_product_dict[a_tup][b_tup]=result
			except KeyError:
				scalar_product_dict[a_tup]={b_tup:result}
			#if a==b and result==0:
			#	print a
			return together(result)
				
					
			
def product_mat_func(i,j):
	return scalar_product(bla_rev[i],bla_rev[j])

def scalar_product_cartan(a,b,index):
	new_vector=highest_weight_vector
	for i in range(len(b)):
		new_vector=subtract(new_vector,vector[b[i]-1])
	return new_vector[index]*scalar_product(a,b)	

def scalar_product_gen(a,b,c):
	first=a[1]
	second=b[1]
	if len(first)!=len(second):
		return 0
	else:
		result=0
		for i in range(0,len(a),2):
			for j in range(0,len(b),2):
				if c=="root":
					result+=a[i]*b[j]*scalar_product(a[i+1],b[j+1])
				else:
					result+=a[i]*b[j]*scalar_product_cartan(a[i+1],b[j+1],c)
		return together(result)
			
		
			

def gram_schmidt(l):
	#print l
	k=[]
	for i in range(len(l)):
		#print i
		h=[1,l[i]]
		for j in k:
			#print "enter"
			coeff=-scalar_product_gen([1,l[i]],j,"root")/scalar_product_gen(j,j,"root")
			for p in range(0,len(j),2):	
				marker=0
				for q in range(0,len(h),2):
					if h[q+1]==j[p+1]:
						h[q]+=coeff*j[p]
						marker=1
				if marker==0:
					h+=[coeff*j[p],j[p+1]]
		k.append(h)
	#print "k is", k
	for i in k:
		#print "this", i
		divide=scalar_product_gen(i,i,"root")
		for j in range(0,len(i),2):
			i[j]=i[j]/sqrt(divide)
	return k				
	
	


def reverse_each(l):
	result=[]
	for i in l:
		new_l=[]
		for j in range(len(i)):
			new_l.append(i[-j-1])
		result.append(new_l)
	return result		

def rep_matrices_func(i,j):
	#print i,j
	ket=states_dict[i+1][0]
	bra=states_dict[j+1][0]
	ket_copy=[]
	for k in range(0,len(ket),2):
		ket_copy+=[ket[k],[now]+ket[k+1]]
	#print ket_copy,bra
	#print i,j,ket_copy,bra
	return scalar_product_gen(bra,ket_copy,"root")

def cartan_gen_func(i,j):
	if i==j:
		ket=states_dict[i+1][0]
		return scalar_product_gen(ket,ket,index)
	else:
		return 0

def simplify_comm(comm):
	if type(comm[1])!=list:
		a=comm[0]
		b=comm[1]
		c=rep_matrices[tuple(0 for j in range(a-1))+(1,)+tuple(0 for j in range(total-a))]
		d=rep_matrices[tuple(0 for j in range(b-1))+(1,)+tuple(0 for j in range(total-b))]
		newmat=c*d-d*c
		for i in range(newmat.shape[0]):
			for j in range(newmat.shape[1]):
				newmat[i*newmat.shape[1]+j]=together(newmat[i*newmat.shape[1]+j])
		return newmat
	else:
		a=comm[0]
		b=simplify_comm(comm[1])
		c=rep_matrices[tuple(0 for j in range(a-1))+(1,)+tuple(0 for j in range(total-a))]
		newmat=c*b-b*c
		for i in range(newmat.shape[0]):
			for j in range(newmat.shape[1]):
				newmat[i*newmat.shape[1]+j]=together(newmat[i*newmat.shape[1]+j])
		return newmat

def f(nlist):
	#print 'nlist', nlist
	if nlist[1]==[]:
	#	print 'trivial'
		return []
	elif type(nlist[1])!=list:
	#	print 'number'
		return nlist
	elif type(nlist[1][1])!=list:
	#	print 'commutator'
		if (nlist[1][0]=='c' and nlist[1][1]=='c'):
			return []
		elif (nlist[1][0]=='c' and nlist[1][1]>0):
			return [prod(nlist[0],sroots[nlist[1][1]-1]),nlist[1][1]]
		
		elif (nlist[1][0]=='c' and nlist[1][1]<0):
			return [prod(-1,prod(nlist[0],sroots[-nlist[1][1]-1])),nlist[1][1]]
		elif (nlist[1][0]>0 and nlist[1][1]=='c'):
			return [prod(-1,prod(nlist[0],sroots[nlist[1][0]-1])),nlist[1][0]]
		
		elif (nlist[1][0]<0 and nlist[1][1]=='c'):
			return [prod(nlist[0],sroots[-nlist[1][0]-1]),nlist[1][0]]
		elif (nlist[1][0]>0 and nlist[1][1]<0):
			if abs(nlist[1][1])==nlist[1][0]:
				return [prod(nlist[0],sroots[nlist[1][0]-1]),'c']
			else:
				return []
		elif (nlist[1][0]<0 and nlist[1][1]>0):
			if abs(nlist[1][0])==nlist[1][1]:
				return [prod(-1,prod(nlist[0],sroots[nlist[1][1]-1])),'c']
			else:
				return []
		#elif nlist[1][0]>nlist[1][1]:
		#	return [prod(-1,nlist[0]),[nlist[1][1],nlist[1][0]]]
		else:
			return nlist
	else:
	#	print 'normal'
		slist=nlist[1]
	#	print slist
		#nest_level=0
		marker=0
		store_list=[]
		#stop_type=None
		while type(slist[1])==list:
			if (slist[0]<0 and slist[1][0]>0):
				marker=1
				#stop_type='pn-pair'
				break
			elif slist[0]=='c':
				marker=1
				#stop_type='cartan'
				break
			else:
				store_list.append(slist[0])
				#nest_level+=1
				slist=slist[1]
	#	print 'store is', store_list
		if marker==0:
			a=f([1,slist])
			if a!=[1,slist] and a!=[]:
	#			print 'interior evaluation'
				newlist=a[1]
				for i in range(len(store_list)):
					newlist=[store_list[-i-1],newlist]
				return f([prod(nlist[0],a[0]),newlist])
			elif a==[]:
				return []
			elif slist[0]==slist[1]:
	#			print 'zero'
				return []
			else:
	#			print 'cannot be evaluated'
				return nlist
		else:
			result=[]
			#if stop_type=='pn-pair':
	#		print 'step 1'
			scomm=f([1,[slist[0],slist[1][0]]])
			if scomm!=[]:
				first=[scomm[1],slist[1][1]]
				first_coeff=prod(scomm[0],nlist[0])
				thing1=f([1,first])
				for i in range(0,len(thing1),2):
					individual1=thing1[i+1]
					storeind1=individual1
					for j in range(len(store_list)):
						storeind1=[store_list[-j-1],storeind1]
					nres=f([prod(first_coeff,thing1[i]),storeind1])
					if result==[]:
						result=nres
					elif nres==[]:
						pass
					else:
						result=newadd(result,nres)
	#		print 'step 2'
			scomm=f([1,[slist[0],slist[1][1]]])
			if scomm!=[]:
				second=[]
				for k in range(0,len(scomm),2):
					if scomm[k+1]=='c' or (type(scomm[k+1])!=list and scomm[k+1]<0):
						second.append(prod(-1,scomm[k]))
						second.append([scomm[k+1],slist[1][0]])
					else:
						second.append(scomm[k])
						second.append([slist[1][0],scomm[k+1]])
				for i in range(0,len(second),2):
					individual2=second[i+1]
					storeind2=individual2
					for j in range(len(store_list)):
						storeind2=[store_list[-j-1],storeind2]
					nres=f([prod(second[i],nlist[0]),storeind2])
					if result==[]:
						result=nres
					elif nres==[]:
						pass
					else:
						result=newadd(result,nres)
			return result

def prod(a,b):
	if type(a)==list and type(b)==list:
		result=0
		for i in range(len(a)):
			result+=a[i]*b[i]
		return result
	elif type(a)!=list and type(b)!=list:
		return a*b
	elif type(a)!=list and type(b)==list:
		result=[]
		for i in range(len(b)):
			result.append(a*b[i])
		return result
	else:
		result=[]
		for i in range(len(a)):
			result.append(b*a[i])
		return result

def newadd(a,b):
	for i in range(0,len(b),2):
		new=b[i+1]
		marker=0
		for j in range(0,len(a),2):
			if a[j+1]==b[i+1]:
				a[j]=diffadd(a[j],b[i])
				marker=1
				break 
		if marker==0:
			a+=[b[i],b[i+1]]
	return a

def diffadd(a,b):
	if type(a)!=list and type(b)!=list:
		return a+b
	elif type(a)==list and type(b)==list:
		res=[]
		for i in range(len(a)):
			res.append(a[i]+b[i])
		return res
			

def form_comm(chain,com_highest):
	factor=com_highest[0]
	rest=com_highest[1]
	for i in range(len(chain)):
		rest=[-chain[-i-1],rest]
	return [com_highest[0],rest]

def cmat_func(i,j):
	return store_cartan[i][j]

def adj(l):
	if type(l[1][1])!=list:
		if -l[1][0]>-l[1][1]:
			return [l[0],[-l[1][1],-l[1][0]]]
		else:
			return [-l[0],[-l[1][0],-l[1][1]]]
	else:
		x=adj([1,l[1][1]])
		return [-l[0]*x[0],[-l[1][0],x[1]]]

def check_sign(x):
	if x>0:
		return 1
	elif x<0:
		return -1
	else:
		return 0

def convert(nest):
	if nest==[]:
		return []
	elif type(nest[1])!=list:
		return nest
	else:
		res=[]
		newnest=nest
		while type(newnest[1])==list:
			res=[newnest[0],]+res
			newnest=newnest[1]
		res=[newnest[1],newnest[0]]+res
		return res

def checkmult(a,b):
	for i in range(len(a)):
		if a[i]!=0:
			if b[i]==0:
				return False
			else:
				multiple=a[i]/b[i]
	if prod(multiple,b)==a:
		return True
	else:
		return False

		


#-------------------------primary functions (ACCESSED only by USER.PY)-----------------------------------------#


def change_notation(family,dimension):
	if type(dimension)!=int:
			raise Exception, "Dimension entered is not an integer"
	if (family=='a' or family=='b' or family=='c' or family=='d' or family=='e' or family=='f' or family=='g'):
		return [family,dimension] 
	elif family=='su':
		family='a'
		dimension-=1
		#print "Note: su(n) is equivalent to a(n-1)\n"
		return [family,dimension]
	elif family=='so':
		if even(dimension):
			family='d'
			dimension=int(dimension/2.)
			#print "Note: so(2n) is equivalent to d(n)\n"
		else:
			family='b'
			dimension=int((dimension-1)/2.)
			#print "Note: so(2n+1) is equivalent to b(n)\n"
		return [family,dimension]
	elif family=='sp':
		if even(dimension):
			family='c'
			dimension=int(dimension/2.)
			#print "Note: sp(2n) is equivalent to c(n)\n"
		else:
			raise Exception, "the sp family can only have even dimension"
		return [family,dimension]
	if dimension<=0:
		raise Exception, "Error, dimension is not positive"
			

def find_sroots(family,dimension):
	global vector
	global total
	vector=[]
	total=dimension
	if family=='a':
			#total=total-1
		for i in range(dimension):
				#print i
				#print sqrt(i), sqrt(3)
			weight1=[0 for j in range(i-1)]+[-sqrt(i)/sqrt(2*(i+1)) for j in range(step(i))]+[1/sqrt(2*j*(j+1)) for j in range(i+1,dimension+1)]
			weight2=[0 for j in range(i)]+[-sqrt(i+1)/sqrt(2*(i+2))]+[1/sqrt(2*j*(j+1)) for j in range(i+2,dimension+1)]
				#print weight1, weight2
			vector.append(subtract(weight1,weight2))
			
	elif family=='b':
		for i in range(dimension-1):
			vector.append([0 for j in range(i)]+[1,-1]+[0 for j in range(dimension-i-2)])
		vector.append([0 for j in range(dimension-1)]+[1])
			
			
	elif family=='c':
		for i in range(dimension-1):
			weight1=[0 for j in range(i-1)]+[-sqrt(i)/sqrt(2*(i+1)) for j in range(step(i))]+[1/sqrt(2*j*(j+1)) for j in range(i+1,dimension)]+[0]
			weight2=[0 for j in range(i)]+[-sqrt(i+1)/sqrt(2*(i+2))]+[1/sqrt(2*j*(j+1)) for j in range(i+2,dimension)]+[0]
			vector.append(subtract(weight1,weight2))
		weight3=[0 for j in range(dimension-2)]+[-sqrt(2*(dimension-1))/sqrt(dimension),0]
		weight4=[0 for j in range(dimension-1)]+[-sqrt(2)/sqrt(dimension)]
		vector.append(subtract(weight3,weight4))	
		
	elif family=='d':
		if dimension>=3:
			for i in range(dimension-1):
				vector.append([0 for j in range(i)]+[1,-1]+[0 for j in range(dimension-i-2)])
			vector.append([0 for j in range(dimension-2)]+[1,1])
		else:
			raise Exception, "The family d is not a pi system for dimension less than 3"
	elif family=='e':
		if dimension==6:
			vector=[[1,-1,0,0,0,0],[0,1,-1,0,0,0],[0,0,1,-1,0,0],[0,0,0,1,-1,0],[-Rational(1,2),-Rational(1,2),-Rational(1,2),-Rational(1,2),Rational(1,2),sqrt(3)/2],[0,0,0,1,1,0]] 
		elif dimension==7:
			vector=[[0, 0, 0, 0, 0, 0, -1, 1], [0, 0, 0, 0, 0, -1, 1, 0], [0, 0, 0, 0, -1, 1, 0, 0], [0, 0, 0, -1, 1, 0, 0, 0], [0, 0, -1, 1, 0, 0, 0, 0], [0, -1, 1, 0, 0, 0, 0, 0],[Rational(1,2),Rational(1,2),Rational(1,2),Rational(1,2),-Rational(1,2),-Rational(1,2),-Rational(1,2),-Rational(1,2)]]
		elif dimension==8:
			vector=[[0, 0, 0, 0, 0, -1, 1, 0], [0, 0, 0, 0, -1, 1, 0, 0], [0, 0, 0, -1, 1, 0, 0, 0], [0, 0, -1, 1, 0, 0, 0, 0], [0, -1, 1, 0, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 0, 0],[Rational(1,2),-Rational(1,2),-Rational(1,2),-Rational(1,2),-Rational(1,2),-Rational(1,2),-Rational(1,2),Rational(1,2)],[1,1,0,0,0,0,0,0]]
		else: 
			raise Exception, "The e family can either be 6,7 or 8 dimensional"
	elif family=='g':
		if dimension==2:
			vector=[[0,1], [sqrt(3)/2,-Rational(3,2)]]
		else:
			raise Exception, "The g family can only be 2 dimensional"
	elif family=='f':
		if dimension==4:
			vector=[[Rational(1,2),-Rational(1,2),-Rational(1,2),-Rational(1,2)],[0,0,0,1],[0,0,1,-1],[0,1,-1,0]]
		else:
			raise Exception, "the f family can only be 4 dimensional"
	else:
		raise Exception, "The family name you entered is not listed in our dictionary. Please check whether you used small letters"
	return vector

def find_proots(family,dimension,sroots):
	global vector
	global total
	global com_dict
	total=dimension
	vector=sroots
	old={}
	current={}
	com_dict={}
	positive_root_dict={}
	cartan_matrix=[[0 for i in range(total)] for j in range(total)]
	
	proots=[]
	for i in range(total):
		old[tuple(0 for j in range(i))+(1,)+tuple(0 for j in range(total-i-1))]=[i+1]+[[0 for j in range(total)] for k in range(2)]+[[],[]]
		proots.append(tuple(0 for j in range(i))+(1,)+tuple(0 for j in range(total-i-1)))

	n=total
	while len(old)>0:
		for element in old:
			for i in range(total):
				check=check1(element,i+1,old[element][1],old[element][0])
				if old[element][0]<=total:
					cartan_matrix[old[element][0]-1][i]=int(round(-check))
				if check>0.9999:
					new_element=increase(element,i+1)
					try:
						current[new_element][1][i]+=1
						current[new_element][2][i]=int(round(check))-1
						current[new_element][3].append(element)
					except KeyError:
						n=n+1
						proots.append(new_element)
						current[new_element]=[n]+[[0 for j in range(i)]+[old[element][1][i]+1]+[0 for j in range(total-i-1)],[0 for j in range(i)]+[int(round(check))-1]+[0 for j in range(total-i-1)]]+[[element],[]]
						predecessor=current[new_element][3][0]
						direction=i+1
						m=(current[new_element][1][i]+current[new_element][2][i])*Rational(1,2)
						j=-m+current[new_element][1][i]
						vector1=vector[i]
						length=sqrt(vector_product(vector1,vector1))
					
					#print direction, length
					#print length
						com_coeff=sqrt(2)/(length*sqrt((m-j+1)*(m+j)))
					#print length, com_coeff
						com_dict[n]=recur([com_coeff,direction,old[element][0]])
					# constructing the algebra
					old[element][4].append(new_element)
			positive_root_dict[element]=[old[element][0], old[element][1], old[element][2],old[element][3],old[element][4]]
			if positive_root_dict[element][0]>total and positive_root_dict[element][4]==[]:
				highest_root=element
				highest_root_serial_no=positive_root_dict[element][0]
		old.clear()
		old.update(current)
		current.clear()
	if not(family=='a' and dimension==1):
		return [proots,cartan_matrix,com_dict,highest_root,highest_root_serial_no]
	else:
		return [proots,cartan_matrix,com_dict]

def find_fweights(dimension,sroots):
	global vector
	global now
	global total
	total=dimension
	vector=sroots
	root_matrix=Matrix(total,total,root_mat_func)
		#root_matrix=array([[vector[i][j] for i in range(total)] for j in range(total)])
	fund_weight=[]
	for i in range(total):
		now=i
		result=Matrix(total,1,result_mat_func)
		ans=root_matrix.LUsolve(result)
		l=[]
		for j in range(total):
			l.append(ans[j])
		fund_weight.append(l) 
	return fund_weight	

def lowest_root_gen(family,dimension,highest_weight):
	changed_details=change_notation(family,dimension)
	family=changed_details[0]
	dimension=changed_details[1]
	sroots=find_sroots(family,dimension)
	info=find_proots(family,dimension,sroots)
	fweights=find_fweights(dimension,sroots)
	positive_root_dict=info[0]
	cartan_matrix=info[1]
	if not(family=='a' and dimension==1):
		highest_root=info[3]
		highest_root_serial_no=info[4]
	else:
		highest_root=None
		highest_root_serial_no=None
	return lowest_root_gen_aux(family,dimension,sroots,cartan_matrix,fweights,highest_root,highest_root_serial_no,highest_weight)


def lowest_root_gen_aux(family,dimension,sroots,cartan_matrix,fund_weight,highest_root,highest_root_serial_no,highest_weight):
	global total
	global vector
	total=dimension
	vector=sroots
	if not(family=='a' and dimension==1):
		lowest_root=tuple(-highest_root[i] for i in range(total)) 
		cartan_lowest=[0 for i in range(total)]
		lines_with_sr=[0 for i in range(total)]
		length=calculate_length(lowest_root)
		for i in range(total):
			newcheck=check1(lowest_root,i+1,[0 for j in range(total)],highest_root_serial_no)
			lines=calculate_lines(lowest_root,i+1)
			cartan_lowest[i]=int(round(-newcheck))
			lines_with_sr[i]=lines
			linear_relation=tuple(lowest_root[i]*vector_product(vector[i],vector[i])/length**2 for i in range(total))
		if highest_weight!=None:
		# highest weight procedure for finding representations
			highest_weight_vector=[0 for i in range(total)]
			for i in range(total):
				highest_weight_vector=add(highest_weight_vector,mult(fund_weight[i],highest_weight[i]))
		else:
			highest_weight_vector=None
			
		result_dict={'lowest_root':lowest_root,'cartan_lowest':cartan_lowest,'lines_with_sr':lines_with_sr,'cartan_matrix':cartan_matrix,'linear_relation':linear_relation,'highest_weight_vector':highest_weight_vector,'total':total}
	#print cartan_matrix	
		return result_dict
	else:
		if highest_weight!=None:
		# highest weight procedure for finding representations
			highest_weight_vector=[0 for i in range(total)]
			for i in range(total):
				highest_weight_vector=add(highest_weight_vector,mult(fund_weight[i],highest_weight[i]))
		else:
			highest_weight_vector=None
		#print {'total':1,'cartan_matrix':[2],'highest_weight_vector':highest_weight_vector}
		return {'total':1,'cartan_matrix':[[2,]],'highest_weight_vector':highest_weight_vector}

def find_weights(family,dimension,highest_weight,ad):
	global total
	global scalar_product_dict
	global highest_weight_vector
	global bla_rev
	global states_dict
	global uroots
	global sroots
	global c
	global product_matrix_cart
	global highest_list
	global store_cartan
	changed_details=change_notation(family,dimension)
	family=changed_details[0]
	dimension=changed_details[1]
	if ad==1 and not(family=='a' and dimension==1):
		sroots=find_sroots(family,dimension)
		info=find_proots(family,dimension,sroots)
		highest_serial=info[4]
		comm_dict=info[2]
		#print comm_dict
		keylist=[keys for keys in comm_dict] 
		for keys in keylist:
			comm_dict[-keys]=adj(comm_dict[keys])
		highest_list=com_dict[highest_serial]
	some_dict=lowest_root_gen(family,dimension,highest_weight)
	total=some_dict['total']
	highest_weight_vector=some_dict['highest_weight_vector']
	cartan_matrix=some_dict['cartan_matrix']

	G=DiGraph()

	prim={}
	mod={}
	weight_dict_high={}
	states_dict={}
	scalar_product_dict={}


	prim[tuple(highest_weight[i] for i in range(total))]=[1]+[[highest_weight[j] for j in range(total)]]+[[],1,[[]]]

	n=1
	check=0
	lowest=tuple(highest_weight[i] for i in range(total))
	n1=1
	states_dict[1]=[[1,[]],tuple(highest_weight[i] for i in range(total)),1]

	G.add_node(tuple(highest_weight[i] for i in range(total)),deg=1)
	G.node[tuple(highest_weight[i] for i in range(total))]['store']=tuple(highest_weight[i] for i in range(total))
	G.node[tuple(highest_weight[i] for i in range(total))]['store_deg']=1
	while len(prim)>0:
		for element in prim:
			#print prim[element][0],element, prim[element][1],prim[element][2],prim[element][3]#,prim[element][4]
			weight_dict_high[element]=[prim[element][0], prim[element][1],prim[element][2],prim[element][3],prim[element][4]]
			for i in range(total):
				if element[i]>0:
					#print element,cartan_matrix[i]
					new_element=subtract(element,cartan_matrix[i])
					check=1
				else:
					if prim[element][1][i]>abs(element[i]):
						new_element=subtract(element,cartan_matrix[i])
						check=1
				if check:
					check=0
					l=[]
					for member in prim[element][4]:
						l.append(member+[i+1])
					
					try:
						if prim[element][1][i]==0:
							mod[new_element][1][i]=element[i]
						else:
							mod[new_element][1][i]=prim[element][1][i]
						mod[new_element][2].append([element,i+1])
						mod[new_element][3]+=prim[element][3]
						mod[new_element][4]+=l
						G.add_edge(element,new_element,root_no=i)
						
					except KeyError:
						n=n+1
						lowest=new_element
						if prim[element][1][i]==0: 
							mod[new_element]=[n]+[[0 for j in range(i)]+[element[i]]+[0 for j in range(total-i-1)]]+[[[element,i+1]]]+[prim[element][3]]+[l]
						else:
							mod[new_element]=[n]+[[0 for j in range(i)]+[prim[element][1][i]]+[0 for j in range(total-i-1)]]+[[[element,i+1]]]+[prim[element][3]]+[l]
						G.add_node(new_element)
						G.node[new_element]['store']=new_element
						G.add_edge(element,new_element,root_no=i)
							
		for new_element in mod:
			if ad==0 or (family=='a' and dimension==1):
				bla=mod[new_element][4]
				bla_rev=reverse_each(bla)
				#print new_element, len(bla_rev)
				l=[]
				bla_mod=[]
				if len(bla_rev)>1:
					product_matrix=Matrix(len(bla_rev),len(bla_rev),product_mat_func)
					#if new_element==(0,0,0,0,0):
						#print product_matrix
						#print mod
					product_matrix_rref=product_matrix.rref()[1]
					for i in product_matrix_rref:
						l.append(bla_rev[i])
						bla_mod.append(bla[i])
					#print new_element,len(l)
					mod[new_element][4]=bla_mod	
						
			
				elif len(bla_rev)==1:
					l.append(bla_rev[0])
				m=0
				for j in gram_schmidt(l):
					n1=n1+1
					m=m+1
					states_dict[n1]=[j,new_element,m]
				#if len(bla_rev)>1:
				#	print states_dict
				
				degeneracy=len(l)
				G.node[new_element]['deg']=degeneracy
				G.node[new_element]['store_deg']=degeneracy
			else:
				if new_element!=tuple(0 for k in range(total)):
					#print new_element
					bla=mod[new_element][4]
					bla_rev=reverse_each(bla)
					for each in bla_rev:
						com=form_comm(each,highest_list)
						#print com
						result=f(com)
						#print each, result[1]
						marker=0
						if result==[]:
							pass
						elif type(result[1])!=list:
							chain=each
							check_zero=scalar_product(chain,chain)
							if check_zero!=0:
								if result[0]<0:
									sign=-1
								else:
									sign=1
								marker=1
						else:
							for keys in comm_dict:
								if comm_dict[keys][1]==result[1]:
									#print keys
									chain=each
									check_zero=scalar_product(chain,chain)
									if check_zero!=0:
										if result[0]<0:
											sign=-1*check_sign(comm_dict[keys][0])
										else:
											sign=1*check_sign(comm_dict[keys][0])
										marker=1
										break
						if marker==1:
							break
					if marker==0:
						checker=bla_rev[0]
						com=form_comm(checker,highest_list)
						#print com
						result=f(com)[1]
						check_nest=convert(result)
						check_nest.sort()
						#print check_nest
						for keys in comm_dict:
							ref=comm_dict[keys][1]
							ref_nest=convert(ref)
							ref_nest.sort()
							if ref_nest==check_nest:
								rightkey=keys
								break
						for each in bla_rev:
							com=form_comm(each,highest_list)
						#print com
							result=f(com)
							if result==[]:
								pass
							else:
								chain=each
								check_zero=scalar_product(chain,chain)
								if check_zero!=0:
									nest=convert(result[1])
									first=form_comm(nest,[1,result[1]])
									#print first
									first_res=f(first)[0]
									#print first_res
									second=form_comm(nest,[1,comm_dict[rightkey][1]])
									second_res=f(second)[0]
									for things in range(len(first_res)):
										if first_res[things]!=0:
											sign=check_sign(first_res[things]*second_res[things])
											break
									break
						sign=check_sign(comm_dict[rightkey][0])*sign
					normalize=scalar_product(chain,chain)
					ans=[sign*(1/sqrt(normalize)),chain]
					#print normalize, ans
					n1=n1+1
					states_dict[n1]=[ans,new_element,1]
					G.node[new_element]['deg']=1
					G.node[new_element]['store_deg']=1
				else:
					bla=mod[new_element][4]
					bla_rev=reverse_each(bla)
					store_cartan=[]
					store_bla_rev=[]
					for i in range(len(bla_rev)):
						com=form_comm(bla_rev[i],highest_list)
						test=f(com)
						marker=0
						for l in range(0,len(test),2):
							#print test[l]
							if test[l]!=list:
								if test[l]!=0:
									marker=1
									break
							else:
								for values in test[l]:
									if values!=0:
										marker=1
										break
						if marker==1:
							#if test[l] not in store_cartan:
							#	store_cartan.append(test[l])
							#	store_bla_rev.append(bla_rev[i])
							is_mult=0
							for things in store_cartan:
								if checkmult(things,test[l]):
									is_mult=1
									break
							if is_mult==0:
								store_cartan.append(test[l])
								store_bla_rev.append(bla_rev[i])
					
									
								
					#print store_cartan
					c=Matrix(total, total, cmat_func)
					#print c
					cinv=c.inv()
					m=0
					start=n1
					for i in range(total):
						l=[]
						for j in range(total):
							l+=[cinv[i*total+j],store_bla_rev[j]]
						norm=sqrt(scalar_product_gen(l,l,'root'))
						for k in range(0,len(l),2):
							if norm==0:
								l[k]=0
							else:
								l[k]*=1/norm
						n1+=1
						m+=1
						states_dict[n1]=[l,new_element,m]
					#for num1 in range(n1,start,-1):
					#	for num2 in range(n1,start,-1):
					#		print num1,num2,scalar_product_gen(states_dict[num1][0],states_dict[num2][0],'root')
					
					G.node[new_element]['deg']=total
					G.node[new_element]['store_deg']=total
					
					
						
				
				
				
		 
		prim.clear()
		prim.update(mod)
		mod.clear()
	return [G,states_dict,scalar_product_dict,highest_weight_vector,n1]

def find_matrices(dimension,sroots,proots,dictionary,states_dictionary,scalar_product_dictionary,highest,dimension_rep,total_pos):
	global vector
	global total
	global com_dict
	global states_dict
	global scalar_product_dict
	global highest_weight_vector
	global now
	global index
	global rep_matrices
	highest_weight_vector=highest
	total=dimension
	vector=sroots
	com_dict=dictionary
	states_dict=states_dictionary
	#print states_dict
	scalar_product_dict=scalar_product_dictionary
	rep_matrices={}

	for i in range(total):
		now=i+1
		#print now
		current_mat=Matrix(dimension_rep,dimension_rep,rep_matrices_func)
		#print current_mat
		#newmat=sqrt(1/together(((current_mat.H*current_mat).trace())))*current_mat
		#print (current_mat.H*current_mat).trace()
		newmat=current_mat
		#print newmat
		rep_matrices[proots[i]]=newmat
		
	for elements in com_dict:
		if elements>total:
			mult_list=com_dict[elements]
			current_mat=mult_list[0]*simplify_comm(mult_list[1])
			#newmat=sqrt(1/together(((current_mat.H*current_mat).trace())))*current_mat
			rep_matrices[proots[elements-1]]=current_mat

	#print total_pos
	newlist=[keys for keys in rep_matrices]
	for keys in newlist:
		rep_matrices[tuple(mult(keys,-1))]=rep_matrices[keys].H
	for i in range(total):
		index=i
		cartan_gen=Matrix(dimension_rep,dimension_rep,cartan_gen_func)
		#rep_matrices[i]=sqrt(1/together(((cartan_gen.H*cartan_gen).trace())))*cartan_gen
		rep_matrices[i]=cartan_gen
	return rep_matrices
	
	



	



