from sympy import *

def f(nlist):
	print 'nlist', nlist
	if nlist[1]==[]:
		print 'trivial'
		return []
	elif type(nlist[1])!=list:
		print 'number'
		return nlist
	elif type(nlist[1][1])!=list:
		print 'commutator'
		if (nlist[1][0]=='c' and nlist[1][1]=='c'):
			return []
		elif (nlist[1][0]=='c' and nlist[1][1]>0):
			return [prod(nlist[0],alphadict[nlist[1][1]]),nlist[1][1]]
		
		elif (nlist[1][0]=='c' and nlist[1][1]<0):
			return [prod(-1,prod(nlist[0],alphadict[-nlist[1][1]])),nlist[1][1]]
		elif (nlist[1][0]>0 and nlist[1][1]=='c'):
			return [prod(-1,prod(nlist[0],alphadict[nlist[1][0]])),nlist[1][0]]
		
		elif (nlist[1][0]<0 and nlist[1][1]=='c'):
			return [prod(nlist[0],alphadict[-nlist[1][0]]),nlist[1][0]]
		elif (nlist[1][0]>0 and nlist[1][1]<0):
			if abs(nlist[1][1])==nlist[1][0]:
				return [prod(nlist[0],alphadict[nlist[1][0]]),'c']
			else:
				return []
		elif (nlist[1][0]<0 and nlist[1][1]>0):
			if abs(nlist[1][0])==nlist[1][1]:
				return [prod(-1,prod(nlist[0],alphadict[nlist[1][1]])),'c']
			else:
				return []
		elif nlist[1][0]>nlist[1][1]:
			return [prod(-1,nlist[0]),[nlist[1][1],nlist[1][0]]]
		else:
			return nlist
	else:
		print 'normal'
		slist=nlist[1]
		print slist
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
		print 'store is', store_list
		if marker==0:
			a=f([1,slist])
			if a!=[1,slist] and a!=[]:
				print 'interior evaluation'
				newlist=a[1]
				for i in range(len(store_list)):
					newlist=[store_list[-i-1],newlist]
				return f([prod(nlist[0],a[0]),newlist])
			elif a==[]:
				return []
			elif slist[0]==slist[1]:
				print 'zero'
				return []
			else:
				print 'cannot be evaluated'
				return nlist
		else:
			result=[]
			#if stop_type=='pn-pair':
			print 'step 1'
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
					elif nres[1]!=result[1]:
						return Exception, 'doesnt work' 
					else:
						result[0]=add(result[0],nres[0])
			print 'step 2'
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
					elif nres[1]!=result[1]:
						return Exception, 'doesnt work' 
					else:
						result[0]=add(result[0],nres[0])
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

def add(a,b):
	if type(a)==list and type(b)==list:
		res=[]
		for i in range(len(a)):
			res.append(a[i]+b[i])
		return res
	elif type(a)!=list and type(b)!=list:
		return a+b
	else:
		return Exception, 'different types added'

alphadict={1:[1,0], 2: [-Rational(1,2),sqrt(3)/2]}

print f([1,[-2,[-1,[-1,[1,2]]]]])

	
			
				
					
							
						
			
