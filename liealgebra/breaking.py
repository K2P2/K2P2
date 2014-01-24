from grp_rep_func import *
#from math import *
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pltpdf
import re
from networkx import *
from sympy.physics.quantum import TensorProduct

#f=open('mat.out')
#g=open('weight.out')


#---------------------------------------------------SECONDARY FUNCTIONS------------------------------------------#

def generate_dynkin(Class,Dim):
	if type(Dim)==int:
		pass
	else:
		raise Exception, 'dimension must be an integer'
	if Class=='a':
		if Dim<1:
			raise Exception, 'invalid dimension'
		G=path_graph(Dim)
		for i in range(Dim-1):
			G[i][i+1]['type']=1
		return G
	elif Class=='d':
		if Dim<3:
			raise Exception,  'invalid dimension'
		G=path_graph(Dim-2)
		for i in range(Dim-3):
			G[i][i+1]['type']=1
		G.add_nodes_from([Dim-2,Dim-1])
		G.add_edges_from([(Dim-3,Dim-2),(Dim-3,Dim-1)])
		G[Dim-3][Dim-2]['type']=1
		G[Dim-3][Dim-1]['type']=1
		return G
	elif Class=='b':
		if Dim==1:
			return generate_dynkin('a',1)
		elif Dim<1:
			raise Exception, 'invalid dimension'
		else:
			G=path_graph(Dim-1)
			for i in range(Dim-2):
				G[i][i+1]['type']=1
			G.add_node(Dim-1)
			G.add_edge(Dim-2,Dim-1)
			G[Dim-2][Dim-1]['type']=2
			for i in range(Dim-1):
				G.node[i]['length']='big'
			G.node[Dim-1]['length']='small'
			return G
	elif Class=='c':
		if Dim==1:
			return generate_dynkin('a',1)
		elif Dim<1:
			raise Exception, 'invalid dimension'
		else:
			G=path_graph(Dim-1)
			for i in range(Dim-2):
				G[i][i+1]['type']=1
			G.add_node(Dim-1)
			G.add_edge(Dim-2,Dim-1)
			G[Dim-2][Dim-1]['type']=2
			for i in range(Dim-1):
				G.node[i]['length']='small'
			G.node[Dim-1]['length']='big'
			return G
	elif Class=='g':
		if Dim==2:
			G=path_graph(2)
			G[0][1]['type']=3
			G.node[0]['length']='small'
			G.node[1]['length']='big'
			return G
		else:
			raise Exception, "The g family can only be 2 dimensional"
	elif Class=='f':
		if Dim==4:
			G=path_graph(4)
			G[0][1]['type']=1
			G[1][2]['type']=2
			G[2][3]['type']=1
			G.node[0]['length']='small'
			G.node[1]['length']='small'
			G.node[2]['length']='big'
			G.node[3]['length']='big'
			return G
		else:
			raise Exception, "the f family can only be 4 dimensional"
	elif Class=='e':
		if Dim==6 or Dim==7 or Dim==8:
			G=path_graph(Dim-1)
			for i in range(Dim-2):
				G[i][i+1]['type']=1
			G.add_node(Dim-1)
			G.add_edge(Dim-4,Dim-1)
			G[Dim-4][Dim-1]['type']=1
			return G
		else:
			raise Exception, "the f family can only be 4 dimensional"
	elif Class=='su':
		if Dim>1:
			Class='a'
			Dim-=1
			#print "Note: su(n) is equivalent to a(n-1)\n"
			return generate_dynkin(Class,Dim)
		else:
			raise Exception, 'invalid dimension'
	elif Class=='so':
		if Dim%2==0 and Dim>=6:
			Class='d'
			Dim=Dim/2
			#print Class,Dim
			return generate_dynkin(Class,Dim)
		elif Dim%2!=0 and Dim>=3:
			Class='b'
			Dim=(Dim-1)/2
			return generate_dynkin(Class,Dim)
		else:
			raise Exception, 'invalid dimension'
	elif Class=='sp':
		if Dim%2==0 and Dim>=2:
			Class='c'
			Dim=Dim/2
			return generate_dynkin(Class,Dim)
		else:
			raise Exception, 'invalid dimension'
	else:
		raise Exception, "invalid class name"


def draw_graph(a):
	e1=[(u,v) for (u,v,d) in a.edges(data=True) if d['type']==1]
	e2=[(u,v) for (u,v,d) in a.edges(data=True) if d['type']==2]
	e3=[(u,v) for (u,v,d) in a.edges(data=True) if d['type']==3]

	pos=nx.spring_layout(a)

	# nodes
	nx.draw_networkx_nodes(a,pos)
	nx.draw_networkx_labels(a,pos)
	
	# edges
	nx.draw_networkx_edges(a,pos,edgelist=e1)
	
	nx.draw_networkx_edges(a,pos,edgelist=e2,alpha=0.5,edge_color='b',style='dashed')
	nx.draw_networkx_edges(a,pos,edgelist=e3,alpha=0.5,edge_color='g',style='dotted')
	#plt.show()

def draw_breaking(a):
	e1=[(u,v) for (u,v,d) in a.edges(data=True) if d['mode']==1]
	e2=[(u,v) for (u,v,d) in a.edges(data=True) if d['mode']==2]

	pos=nx.spring_layout(a)

	# nodes
	nx.draw_networkx_nodes(a,pos,node_shape='s',node_size=700,node_color='b',alpha=0.5,linewidths=None)
	nx.draw_networkx_labels(a,pos)
	
	# edges
	nx.draw_networkx_edges(a,pos,edgelist=e1)
	
	nx.draw_networkx_edges(a,pos,edgelist=e2,alpha=0.5,edge_color='b')
	#nx.draw_networkx_edges(a,pos,edgelist=e3,alpha=0.5,edge_color='g',style='dotted')
	#plt.show()

def identify_dynkin(G):
	l=len(G)
	#print l
	em=isomorphism.categorical_edge_match('type',1)
	nm=isomorphism.categorical_node_match('length','small')
	H=generate_dynkin('a',l)
	GM=isomorphism.GraphMatcher(H,G,edge_match=em)
	#print H.nodes(),G.nodes(),H.edges(),G.edges()
	if GM.is_isomorphic():
		#print True
		return {'group':['su',l+1],'mapping':GM.mapping}
	if l>=2:
		H=generate_dynkin('b',l)
		GM=isomorphism.GraphMatcher(H,G,node_match=nm,edge_match=em)
		if GM.is_isomorphic():
			return {'group':['so',2*l+1],'mapping':GM.mapping}
		H=generate_dynkin('c',l)
		GM=isomorphism.GraphMatcher(H,G,node_match=nm,edge_match=em)
		if GM.is_isomorphic():
			return {'group':['sp',2*l],'mapping':GM.mapping}
		if l==2:
			H=generate_dynkin('g',2)
			GM=isomorphism.GraphMatcher(H,G,node_match=nm,edge_match=em)
			if GM.is_isomorphic():
				return {'group':['g',2],'mapping':GM.mapping}
		if l>=3:
			H=generate_dynkin('d',l)
			GM=isomorphism.GraphMatcher(H,G,edge_match=em)
			if GM.is_isomorphic():
				return {'group':['so',2*l],'mapping':GM.mapping}
			if l==4:
				H=generate_dynkin('f',4)
				GM=isomorphism.GraphMatcher(H,G,node_match=nm,edge_match=em)
				if GM.is_isomorphic():
					return {'group':['f',4],'mapping':GM.mapping}
			if l==6:
				H=generate_dynkin('e',6)
				GM=isomorphism.GraphMatcher(H,G,edge_match=em)
				if GM.is_isomorphic():
					return {'group':['e',6],'mapping':GM.mapping}
			if l==7:
				H=generate_dynkin('e',7)
				GM=isomorphism.GraphMatcher(H,G,edge_match=em)
				if GM.is_isomorphic():
					return {'group':['e',7],'mapping':GM.mapping}
			if l==8:
				H=generate_dynkin('e',8)
				GM=isomorphism.GraphMatcher(H,G,edge_match=em)
				if GM.is_isomorphic():
					return {'group':['e',7],'mapping':GM.mapping}
		
	raise Exception, "invalid diagram"

def remove_edges(graph,rem):
	for i in graph.edges(data=True):
		if i[2]['root_no']==rem:
		        graph.remove_edge(i[0],i[1])
	return graph

def change_node(graph,number,rem):
	new_map={number:rem}
	graph1=relabel_nodes(graph,new_map)
	return graph1

def change_mapping(mapping,number,rem):
	for entry in mapping:
		if mapping[entry]==number:
			#print 'yes'
			mapping[entry]=rem
	return mapping

def change_node_tuples(graph,rem,relation,mapping,total):
	newlist=[0 for i in range(total)]
	for i in range(len(relation)):
		newlist[mapping[i]]=relation[i]
	#print newlist
	newdict={}
	#print 'rem is', rem
	#if (0,0,0,0,-1,0) in graph:
		#print 'yes'
	for nodes in graph:
		newvalue=0
		for i in range(total):
			newvalue+=newlist[i]*nodes[i]
		newdict[nodes]=tuple(nodes[i] for i in range(rem))+(newvalue,)+tuple(nodes[i] for i in range(rem+1,total,1))
	for i in newdict:
		if newdict[i]==(0,0,1,0,-1,0):
			#print "pred is", i
			break
	graph1=relabel_nodes(graph,newdict)
	return graph1

def match_spec(nodes,new,newlist):
	logical=True
	for i in range(len(nodes)):
		if i in newlist:
			if nodes[i]!=new[i]:
				logical=False
	return logical

def find_tuple(graph,new,mapping):
	newlist=[]
	for keys in mapping:
		newlist.append(mapping[keys])
	anotherlist=[]
	for nodes in graph:
		if match_spec(nodes,new,newlist):
			anotherlist.append(nodes)
	return anotherlist
			

def join_with_new(graph,rem,cartan,mapping,total):
	#print cartan
	newcartan=[0 for i in range(total)]
	for i in range(len(cartan)):
		newcartan[mapping[i]]=cartan[i]
	newcartan[rem]=2
	#newcartan[0]=-1
	for nodes in graph:
		if nodes[rem]>0:
			start=nodes
			for i in range(nodes[rem]):
				new=subtract(start,newcartan)
				#print new,start
				if new in graph:
					graph.add_edge(start,new,root_no=rem)
				else:
					#raise Exception, "New node not found. Check the algorithm!"
					#pass
					newlist=find_tuple(graph,new,mapping)
					if len(newlist)>1:
						raise Exception, "multiple connections while going down. Please modify the code slightly"
					else:
						new=newlist[0]
						graph.add_edge(start,new,root_no=rem)
				start=new
	return graph


					
	
def first_highest(graph):
	nei_list=[]
	for nodes in graph:
		nei_list+=graph.neighbors(nodes)
	#print nei_list
	for nodes in graph:
		if nodes not in nei_list:
			result=nodes
			break
	try:
		return result
	except	UnboundLocalError:
		return None


def get_nodes(dynlist):	
	if dynlist==[]:
		return [()]
	else:
		l=[]
		m=dynlist[1:]
		for k in get_nodes(m):
			for i in dynlist[0]:	
				l.append(((i,))+k)
		return l

def form_graph(dynlist,mapping):
	l=[() for k in range(len(dynlist))]
	#numlist=[len(dynlist[k]) for k in range(len(dynlist))]
	for i in range(len(dynlist)):
		for nodes in dynlist[i]:
			l[i]+=((nodes,))
	#print l
	newlist=get_nodes(l)
	newgraph=DiGraph()
	newgraph.add_nodes_from(newlist)
	#print newgraph.nodes()
	for nodes in newgraph:
		#rank=[]
		#for i in range(len(nodes)):
		#	for j in range(len(l[i])):
		#		if l[i][j]==nodes[i]:
		#			rank.append(j)
		#			break
		#newgraph.node[nodes]['rank']=rank
		for i in range(len(nodes)):
			nei=dynlist[i].neighbors(nodes[i])
			for j in nei:
				newnode=tuple(nodes[k] for k in range(i))+(j,)+tuple(nodes[k] for k in range(i+1,len(nodes),1))
				num=mapping[i][dynlist[i].edge[nodes[i]][j]['root_no']]
				newgraph.add_edge(nodes,newnode,root_no=num)
	for nodes in newgraph:
		deg=1
		for i in range(len(nodes)):
			deg*=dynlist[i].node[nodes[i]]['deg']
		newgraph.node[nodes]['deg']=deg
	#print newgraph.nodes()
	#print newgraph.edges()
	return newgraph

def form_trial_mat(vec_list,trial):
	l1=[]
	for vec in vec_list:
		l2=[]
		for num in range(len(vec)):
			if num!=trial:
				l2.append(vec[num])
		l1.append(l2)
	return Matrix(l1)

def form_trial_vec(vec_list,trial):
	l1=[]
	for i in range(len(vec_list)):
		l1.append([-vec_list[i][trial]])
	return Matrix(l1)

def form_col_vec(nodes,states_dict):
	l=[]
	for i in range(len(states_dict)):
		if states_dict[i+1][1]==nodes:
			l.append([1])
		else:
			l.append([0])
	return Matrix(l)


def tuple_to_col_vec(tup):
	l=[]
	for i in tup:
		l.append([i])
	return Matrix(l)

def normalize(col_mat):
	a=col_mat.T*col_mat
	return (1/sqrt(together(a[0])))*col_mat

def gram_schmidt_orthogonal(newstate,statelist):
	#print "the bunch is", statelist
	ans=newstate
	#print "the newstate is", newstate
	for states in statelist:
		#print newstate.T,states
		a=newstate.T*states
		#print a[0],together(a[0])
		#print "a is", a[0]
		ans=ans-together(a[0])*states
		#for i in range(newans.shape[0]):
		#	newans[i]=together(newans[i])
		#ans=newans
	#print "the final state is", ans
	prod=ans.T*ans
	if abs(prod[0])<=1.e-7:
		return None
	else:
		#for i in range(ans.shape[0]):
		#	ans[i]=together(ans[i],5)
		return ans

def form_col_vec_ind(states_dict,indices,i):
	l=[]
	for j in range(len(states_dict)):
		if j==indices[i]:
			l.append([1])
		else:
			l.append([0])
	return Matrix(l)

def gram_schmidt_orthogonalize_many(first,second):
	l=[]
	start=first
	for vec in second:
		newvec=vec
		#print "newvec is", newvec
		for other_vec in start:
			a=other_vec.T*vec
			newvec=newvec-together(a[0])*other_vec
		#for i in range(newvec.shape[0]):
			#newvec[i]=together(newvec[i],5)
			
		start.append(newvec)
		#print "start is", start
		l.append(newvec)
	#print first,second,l
	return l

def tens_prod_higgs(higgs,dimlist,index):
	l=[]
	for i in range(len(dimlist)):
		if i==index:
			l.append(higgs)
		else:
			m=[[1]]
			for j in range(1,dimlist[i],1):
				m.append([0])
			l.append(Matrix(m))
	ans=l[0]
	for i in range(1,len(l),1):
		ans=TensorProduct(ans,l[i])
		for j in range(ans.shape[0]):
			for k in range(ans.shape[1]):
				ans[j*ans.shape[1]+k]=together(ans[j*ans.shape[1]+k])
	return ans

def tens_prod(mat,dimlist,index):
	l=[]
	for i in range(len(dimlist)):
		if i==index:
			l.append(mat)
		else:
			l.append(eye(dimlist[i]))
	ans=l[0]
	for i in range(1,len(l),1):
		ans=TensorProduct(ans,l[i])
		for j in range(ans.shape[0]):
			for k in range(ans.shape[1]):
				ans[j*ans.shape[1]+k]=together(ans[j*ans.shape[1]+k])
	return ans
		

#------------------------------------------------PRIMARYFUNCTIONS-------------------------------------------------#

def break_algebra(fam,dim):
	total=dim
	first=generate_dynkin(fam,dim)
	for nodes in first:
		first.node[nodes]['linear_comb']=tuple(0 for i in range(nodes))+(1,)+tuple(0 for i in range(total-nodes-1))
	process={0:[first]}
	info=identify_dynkin(first)['group']
	name=info[0]+'('+str(info[1])+')'
	process_name={0:[name]}
	process_info={0:[info]}
	first.graph['name']=name
	break_info_parent={}
	break_info_child={}

	G=Graph()
	G.add_node(first,height=0,code=info)

	#print G.node[first]['height']

	count=0
	number=-1
	network_info=[[0,None]]

	while True:
		print "Step %s" % str(count)
		print "You can break any of the following algebra/subalgebras in this step"
		print process_name[count]
		print "To break an algebra/subalgebra, enter the position index of the algebra/subalgebra in the above list, enter e to exit"
		action=raw_input('Please enter your choice:')
		if action=='e':
			break
		else:
			try:
				if int(action)<=len(process_name[count]):
					i=int(action)
					m,q,l=[],[],[]
					info=process_name[count][i]
					keep_intact=1
					print "This is a %s algebra\n" % info
					while True:
						if fam=='a' and dim==1:
							print "The following actions are allowed:\n[1]Break via removing a circle\n[3]Show the dynkin diagram\n[4]Keep Intact"
						else:	
							print "The following actions are allowed:\n[1]Break via removing a circle\n[2]Break via adding the lowest weight followed by removing a circle\n[3]Show the dynkin diagram\n[4]Keep Intact"
						choice=input("Please enter your choice:")		
						if choice==1:
							rem=input("remove circle no:")
							try:
								H=process[count][i].copy()
								H.remove_node(rem)
								#print "the algebra was broken down to:"
								keep_intact=0
								conn=connected_components(H)
								for j in conn:
									new=H.subgraph(j).copy()
									newinfo=identify_dynkin(new)['group']
									newname=newinfo[0]+'('+str(newinfo[1])+')'
									m.append(newname)
									q.append(newinfo)
									#print newname
									new.graph['name']=newname
									l.append(new)
									G.add_node(new,height=G.node[process[count][i]]['height']+1,code=newinfo)	
									G.add_edge(process[count][i],new,mode=1,delete=rem)
						#network=remove_edges(network,rem)
						#print network.edges(data=True)
								network_info.append([1,rem])
							except NetworkXError:
								print "problem removing node, please take a look at the dynkin diagram (option [3]) to make sure that the node specified for removing actually exists"
							break
						elif choice==2:
							if (fam=='a' and dim==1):
								print "Error: Invalid input"
							else:
								P=re.compile('\(')
								M=P.search(info)
								family=info[:M.start()]
								P=re.compile('[0-9]+')
								M=P.search(info)
								dimension=int(M.group())
								lowest_dict=lowest_root_gen(family,dimension,None)
								H=process[count][i].copy()
								H.add_node(number)
								lines_with_sr=lowest_dict['lines_with_sr']
								mapping=identify_dynkin(process[count][i])['mapping']
								for k in range(len(lines_with_sr)):
									if lines_with_sr[k]!=0:
										H.add_edge(mapping[k],number)
										H[mapping[k]][number]['type']=lines_with_sr[k]	
								print "Showing the resulting dynkin diagram now. You will be asked to remove a circle subsequently"
								linear_comb=lowest_dict['linear_relation']
								newvector=tuple(0 for k in range(total))
								for j in range(len(linear_comb)):
									newvector=add(newvector,tuple(mult(H.node[mapping[j]]['linear_comb'],linear_comb[j])))
								H.node[number]['linear_comb']=newvector	
								draw_graph(H)
								plt.show()
								rem=input("remove circle no:")
								H.remove_node(rem)
								H=change_node(H,number,rem)
								mapping=change_mapping(mapping,number,rem)
								number-=1
							#print "the algebra was broken down to:"
								keep_intact=0
								conn=connected_components(H)
								for j in conn:
									new=H.subgraph(j).copy()
									newinfo=identify_dynkin(new)['group']
									newname=newinfo[0]+'('+str(newinfo[1])+')'
									m.append(newname)
									q.append(newinfo)
									#print newname
									new.graph['name']=newname
									#print 'new name is', new.graph['name']
									l.append(new)
									#print process[count][i]
									#print count,i
									G.add_node(new,height=G.node[process[count][i]]['height']+1,code=newinfo)		
									G.add_edge(process[count][i],new,mode=2,delete=rem)
							#network=remove_edges(network,rem)
							#network=change_node_tuples(network,rem,lowest_dict['linear_relation'],mapping,total)
							#commutator=lowest_dict['commutator']
							#current_mat=commutator[0]*simplify_comm_new(commutator[1],main['matrices'],mapping)	
							#current_mat=current_mat.H
							#new_mat=sqrt(1/((current_mat.H*current_mat).trace()))*current_mat
							#network=join_with_new(network,rem,lowest_dict['cartan_lowest'],mapping,total)
							#main['matrices'][rem]=new_mat
								network_info.append([2,[rem,lowest_dict['linear_relation'],mapping,total,lowest_dict['cartan_lowest']]])
								break
						elif choice==3:
							print "The dynkin diagram will now be displayed. Please close the graph window to proceed"
							draw_graph(process[count][i])
							plt.show()
						elif choice==4:
							break
						else:
							print "Error: Invalid Input"
					if keep_intact==1:
						pass
					else:
						break_info_parent[count]=[0 for j in range(i)]+[1,]+[0 for j in range(len(process[count])-i-1)]
						break_info_child[count+1]=[0 for j in range(i)]+[1 for j in range(len(l))]+[0 for j in range(len(process[count])-i-1)]
						process[count+1]=process[count][:i]+l+process[count][i+1:]
						process_name[count+1]=process_name[count][:i]+m+process_name[count][i+1:]
						process_info[count+1]=process_info[count][:i]+q+process_info[count][i+1:]
						count+=1
			except ValueError:
				print "invalid input, try again"		
			
			
	return [count,process,process_info,break_info_parent,break_info_child,network_info]

def break_rep(family,dimension,highest_weight,count,process,process_info,break_info_parent,break_info_children,network_info,ad):

	# global declarations 
	#	1. total is the root space dimension of the algebra under consideration
	#	2. trial appears in computation of U(1) elements. For full explanation see that part of the code.
	#	3. vec_list appears in the computation of U(1) elements. It contains the simple roots (expressed as vectors in root space) corresponding to the remaining circles (after removal of a circle via ROUTE 1 of breaking) . (This is global because I use a matrix generating function outside this routine which uses vec_list as input) 

	global vec_list	
	global trial
	global total


	# get all algebra related information 
	#	1. root space dimension
	#	2. simple roots as vectors
	#	3. positive roots -> dynkin labels for non-cartan generators
	#	4. commutation dictionary -> commutation relation between positive roots
	#	5. states dictionary -> states of the representation expressed as lowering operators acting on the highest weight state
	#	6. scalar product dictionary -> dictionary of scalar products between all (including un-orthogonalized) states that appear in the comutation for states
	# 	7. highest weight as a vector
	#	8. dimension of the representation 
	# 	9. matrices for the generators of the algebra, keys are numbers for cartan generators and dynkin labels for others
	#	10. network is the graph of states (weight diagram)

	total=dimension
	sroots=find_sroots(family,dimension)
	output=find_weights(family,dimension,highest_weight,ad)
	info=find_proots(family,dimension,sroots)
	proots=info[0]
	com_dict=info[2]
	pdim=len(info[0])
	states_dict=output[1]
	scalar_product_dict=output[2]
	highest_weight_vector=output[3]
	dimension_rep=output[4]
	mat=find_matrices(dimension,sroots,proots,com_dict,states_dict,scalar_product_dict,highest_weight_vector,dimension_rep,pdim)
	network=output[0]
	network_store=network.copy()

	# matlist is a list of keys and generator matrices -> same info as mat (why do I need this?)

	matlist=[]
	for keys in mat:
		matlist.append([keys,mat[keys]])

	# adim -> adjoint dimension which is the total number of generators 
	adim=len(mat)

	# initialize all the data to be returned
	#	1. particle content will contain info about the broken representations and the states contained in them(as linear combination of original states)
	#	2. U1_vectors is ?
	#	3. U1_gen is a dictionary that returns the U1 generator created when route 1 (simply remove a circle) is employed for breaking 
	#	4. gen_dict is a dicitonary of other generators (details need to be spelled out)
	#	5. higgs_info is information of singlets left out after breaking 

	particle_content={}
	U1_vectors={}
	U1_gen={}	
	gen_dict={}
	higgs_info={}

	# What are these?

	correspond={}
	statelist=[]

	# this has something to do with degenerate nodes (but what exactly?)
	degnodelist=[]
	

	# START OF THE ROUTINE PROPER

	for index in range(count+1):			# Iteration over each breaking step

		# EACH BREAKING STEP
		print "calculating for step", index
		inf=network_info[index]			# contains info on A. Route of breaking B. Node removed


		# PHASE 1



		# IF ROUTE 1 (SIMPLY REMOVE A CIRCLE) IS CHOSEN

		if inf[0]==1:

			# STEP 1 : REMOVE THE CONNECTIONS CORRESPONDING TO DELETED CIRCLE IN THE WEIGHT DIAGRAM
	
			network=remove_edges(network,inf[1])	


			# STEP 2 : FIND THE U(1) ELEMENT  -> Find a vector V orthogonal to simple roots (expressed as vectors in the root space) corresponding to remaining circles. V.H is the U(1) generator as it commutes with all the remaining simple roots (simple raising and lowering generators).

			# storing the vectors (in root space) corresponding to remaining circles in vec_list
			vec_list=[]
			for dia in process[index]:
				for nodes in dia:
					new_vec=[0 for dim in range(total)]
					for dim in range(total):
						new_vec=add(new_vec,mult(sroots[dim],dia.node[nodes]['linear_comb'][dim])) # some of the circles might be linear combination of the initially chosen simple roots. this happens if ROUTE 2 has been employed earlier in the breaking process. In this case the lowest root that is added is a linear combination of initially chosen simple roots.
					vec_list.append(new_vec)

			# The newly chosen U(1) element must also commute with the U(1) elements chosen in earlier steps. The vectors V of earlier steps are stored in U1_vectors. This is now added to vec_list.
			for keys in U1_vectors:
				vec_list.append(U1_vectors[keys])

			#setting up the equation for finding the vector V orthogonal to every vector in vec_list (While doing this we try all solutions which have 1 as one of the components. this is done to ignore overall scaling which leads to infinite number of solutions. however it is not known a priori whether a certain component is 0 or non-zero. we can only set a non zero component to 1. this is why all components are tried one by one. this is denoted by the variable trial. in case we have chosen a non zero component, the process yeilds a solution. otherwise it yeilds an error. 

			sol=0
			for trial in range(total):
				A=form_trial_mat(vec_list,trial)
				C=form_trial_vec(vec_list,trial)
				try:
					B=A.LUsolve(C)
					B1=[B[j] for j in range(trial)]+[1]+[B[j] for j in range(trial,total-1,1)]
					sol=1
					length=0
					for num in B1:
						length+=num**2
					B1=mult(B1,1/sqrt(length)) # this is V normalized

					# V IS FOUND

					U1_vectors[index]=B1
					newmat=zeros(dimension_rep)
					for num in range(len(B1)):
						newmat+=B1[num]*mat[num]
					for i in range(newmat.shape[0]):
						for j in range(newmat.shape[1]):
							newmat[i*newmat.shape[1]+j]=together(newmat[i*newmat.shape[1]+j])

					# U(1) ELEMENT IS FOUND (by simply computing V.H where info on H comes from mat) 

					U1_gen[index]=newmat
				except (ValueError or ZeroDivisionError):
					pass 
				if sol==1:
					break
			if sol==0:
				return Exception, "No U(1) elements found"   # this is a fatal error!

		# IF ROUTE 2 (LOWEST ROOT ADDED AND A CIRCLE IS REMOVED) IS CHOSEN

		elif inf[0]==2:

			# STEP 1 : REMOVE THE CONNECTIONS CORRESPONDING TO DELETED CIRCLE IN THE WEIGHT DIAGRAM		

			network=remove_edges(network,inf[1][0])

			# STEP 2: CHANGE THE DYNKIN LABELS OF THE STATES IN THE WEIGHT DIAGRAM. THE PLACE CORRESPONDING TO THE DELETED CIRCLE IS TAKEN BY THE LOWEST ROOT THAT HAS BEEN ADDED
			network=change_node_tuples(network,inf[1][0],inf[1][1],inf[1][2],inf[1][3])

			# STEP 3 : ADDITIONAL CONNECTIONS CORRESPONDING TO THE LOWEST ROOTS ARE MADE IN THE WEIGHT DIAGRAM

			network=join_with_new(network,inf[1][0],inf[1][4],inf[1][2],inf[1][3])
		else:
			pass

		# PHASE 2
		
		graphlist=process[index]
		infolist=process_info[index]
		maplist=[]
		particle_content_list=[]
		newnetwork=network.copy()
		newnetwork_store=newnetwork.copy()
		gen_dict[index]=[]
		for i in range(len(graphlist)):
			newmap=identify_dynkin(graphlist[i])['mapping']
			maplist.append(newmap)
			newnames=change_notation(process_info[index][i][0],process_info[index][i][1])
			newfam,newdim=newnames[0],newnames[1]
			newsroots=find_sroots(newfam,newdim)
			newproots=find_proots(newfam,newdim,newsroots)[0]
			l=[]
			for things in newproots:
				#print newproots
				newtup=tuple(0 for k in range(total))
				for j in range(len(things)):
					#print graphlist[i].node[newmap[j]]['linear_comb'],things[j]
					newtup=add(newtup,tuple(mult(graphlist[i].node[newmap[j]]['linear_comb'],things[j])))
				l.append(newtup)
			#print l
			list_to_gs=[tuple_to_col_vec(graphlist[i].node[points]['linear_comb']) for points in graphlist[i]]
			newlist=gram_schmidt_orthogonalize_many([list_to_gs[0]],list_to_gs[1:])
			m=[]			
			for things in newlist:
				things=normalize(things)
				newmat=zeros(dimension_rep)
				for j in range(things.shape[0]):
					newmat+=things[j]*mat[j]
				m.append(newmat)
			mat_dict={}
			for things in l:
				for objects in matlist:
					if objects[0]==things:
						mat_dict[things]=objects[1]
						mat_dict[tuple(mult(things,-1))]=objects[1].H
						break
			p=0
			for things in m:
				mat_dict[p]=things
				p+=1
			gen_dict[index].append(mat_dict)
		#finding the generators corresponding to an algebra
			
		while len(newnetwork)>0:
			a=first_highest(newnetwork)
			#print newnetwork_store.node[a]['store']
			#print a
			if a==None:
				raise Exception, "no highest weights found during decomposition"	
			highest_list=[]
			for i in maplist:
				if len(i)==1:
					sublist=[a[i[0]],]
					#if a[i[0]]!=0:
					#	singlet=0
				else:
					sublist=[]
					for nodes in i:
						sublist.append(a[i[nodes]])
						#if a[i[nodes]]!=0:
						#	singlet=0
				highest_list.append(sublist)
			dynlist=[]
			#print highest_list
			for i in range(len(infolist)):
				#print infolist
				new=change_notation(infolist[i][0],infolist[i][1])
				family=new[0]
				dim=new[1]
				highest=highest_list[i]
				print family,dim
				sr=find_sroots(family,dim)
				print sr
				hr=find_proots(family,dim,sr)[0][-1]
				res=tuple(0 for i in range(dim))	
				for i in range(len(hr)):
					res=add(res,mult(sr[i],hr[i]))
				ad=[0 for i in range(dim)]
				for i in range(dim):
					ad[i]=2*vector_product(res,sr[i])/vector_product(sr[i],sr[i])
				if highest==ad:
				#print family,dim,highest
					dynlist.append(find_weights(family,dim,highest,1)[0])
				else:
					dynlist.append(find_weights(family,dim,highest,0)[0])
			newgraph=form_graph(dynlist,maplist)
			#plt.figure(30)
			#draw(newgraph)
			#plt.savefig('test.png')
			#break
			highest_tuple=()
			for i in highest_list:
				highest_tuple+=(tuple(i[k] for k in range(len(i))),)
			#node_dict={a:highest_tuple}
			#print newgraph.edges(data=True)
			#print maplist
			em=isomorphism.categorical_edge_match('root_no',-1)
			GM=isomorphism.DiGraphMatcher(newnetwork,newgraph,edge_match=em)
			logical=False
			for i in GM.subgraph_isomorphisms_iter():
				if a in i and i[a]==highest_tuple:
					newmap=i
					logical=True
					#if singlet==1:
					#	for keys in newmap:
					#		print newnetwork.node[keys]['store']
					break
			#draw(newgraph)
			#plt.show()
			if logical==False:
				raise Exception, "no isomorphisms where found. possible error in algorithm. please try debugging"	
			#print newma
			revmap={}
			for keys in newmap:
				revmap[newmap[keys]]=keys		
			nodes_matched=[i for i in newmap]
			start=[newmap[a]]
			if newgraph.neighbors(newmap[a])!=[]:
				if newnetwork_store.node[a]['deg']>1:
					try:
						veclist=newnetwork_store.node[a]['vectors']
						act_node=newnetwork_store.node[a]['store']
						span=newnetwork_store.node[a]['deg']
						indices=[]
						for states in states_dict:
							if states_dict[states][1]==act_node:
								indices.append(states-1)	
						for i in range(span):
							newvec=form_col_vec_ind(states_dict,indices,i)
							ans=gram_schmidt_orthogonal(newvec,newnetwork_store.node[a]['vectors'])
							if ans!=None:
								vec=ans
								break
						#print vec
						vec=normalize(vec)
						particle_content[index].append([highest_list,[vec]])
						newnetwork_store.node[a]['vectors'].append(vec)
					except KeyError:
						try:
							particle_content[index].append([highest_list,[form_col_vec(newnetwork_store.node[a]['store'],states_dict)]])
							newnetwork_store.node[a]['vectors']=[form_col_vec(newnetwork_store.node[a]['store'],states_dict)]
						except KeyError:
							particle_content[index]=[[highest_list,[form_col_vec(newnetwork_store.node[a]['store'],states_dict)]]]
							newnetwork_store.node[a]['vectors']=[form_col_vec(newnetwork_store.node[a]['store'],states_dict)]
				
			

				else:					
					try:
						particle_content[index].append([highest_list,[form_col_vec(newnetwork_store.node[a]['store'],states_dict)]])
						newnetwork_store.node[a]['vectors']=[form_col_vec(newnetwork_store.node[a]['store'],states_dict)]
					except KeyError:
						particle_content[index]=[[highest_list,[form_col_vec(newnetwork_store.node[a]['store'],states_dict)]]]
						newnetwork_store.node[a]['vectors']=[form_col_vec(newnetwork_store.node[a]['store'],states_dict)]
			print "creating state dictionary", highest_list
			while start!=[]:
				end=[]
				for things in start:
					for elements in newgraph.neighbors(things):
						if elements not in end:
							end.append(elements)
						#print elements
						if newnetwork_store.node[revmap[elements]]['deg']>1:
							#print 'yes'
							act_element=newnetwork_store.node[revmap[elements]]['store']
							#if act_element==(0,0,-1,1,1):
							#	print "yes"
							try:
								#print newnetwork_store.node[revmap[things]]['vectors']
								#act_parent_list=newnetwork_store.node[revmap[things]]['vectors']
								act_parent_list=[]
								for i in newnetwork_store.node[revmap[things]]['vectors']:
									if i in particle_content[index][-1][1]:
										act_parent_list.append(i)
							except KeyError:
								matvec=form_col_vec(newnetwork_store.node[revmap[things]]['store'],states_dict)
								act_parent_list=[matvec,]
								#print "no_vector"
							for act_parent in act_parent_list:
								#print act_parent_list[0]
								root=newgraph[things][elements]['root_no']
								for dia in graphlist:
									if root in dia:
										combination=dia.node[root]['linear_comb']	
										break
								#print tuple(mult(combination,-1))
								#print "combination is", combination 
								#pos=1
								#for i in range(len(combination)):
								#	if combination[i]<0:
								#		pos=0
								#		break
								#if pos==1:
								#	matrix=mat[combination][1]
								#else:
								#	matrix=mat[combination][1]
								newcomb=tuple(mult(combination,-1))
								for objects in matlist:
									if objects[0]==newcomb:
										matrix=objects[1]
									
								
								#try:
								#	matrix=mat[combination][1]
								#except KeyError:
								#	print combination
								#	newcomb=()
								#	for i in range(len(combination)):
								#			newcomb+=(abs(combination[i]),)
								#matrix=mat[newcomb][0]
								#col_vec=form_col_vec(act_parent,states_dict)
								#print act_parent,matrix,col_vec
								newstate=matrix*act_parent
								#print newstate,matrix,act_parent
								#print revmap[elements],act_parent,newstate
								try:
									#print len(newnetwork_store.node[revmap[elements]]['vectors'])
									#if highest_list==[[0,0,0],[1],[1]]:
									#	print "newstate is", newstate
									#	print "orth is", newnetwork_store.node[revmap[elements]]['vectors']		
									#if highest_list==[[0,0],[1],[1]]:
									#	print "newstate_formed"
									newstate=gram_schmidt_orthogonal(newstate,newnetwork_store.node[revmap[elements]]['vectors'])
									
									if newstate!=None:
									#	if highest_list==[[0,0,0],[1],[1]]:
									#		print "othogonalized is", newstate
										newstate=normalize(newstate)
										#if act_element==(0,0,-1,1,1):
											#print "in"
										newnetwork_store.node[revmap[elements]]['vectors'].append(newstate)
										#if highest_list==[[0,0],[1],[1]]:
											#print "newstate_added"
									#else:
										#if highest_list==[[0,0],[1],[1]]:
										#	print "newstate_exists"
								except KeyError:
									#if act_element==(0,0,-1,1,1):
										#print "in"
									newstate=normalize(newstate)
									newnetwork_store.node[revmap[elements]]['vectors']=[newstate]
								if (newstate not in particle_content[index][-1][1]) and (newstate!=None):
									particle_content[index][-1][1].append(newstate)
						else:
							#if highest_list==[[0,0],[1],[1]]:
								#print "non_degenerate newstate added"
							act_element=newnetwork_store.node[revmap[elements]]['store']
							#if act_element==(0,0,-1,1,1):
								#print "in"
							col_vec=form_col_vec(act_element,states_dict)
							newnetwork_store.node[revmap[elements]]['vectors']=[col_vec]
							if col_vec not in particle_content[index][-1][1]:
								particle_content[index][-1][1].append(col_vec)
				#print highest_list,start,end
				#if highest_list==[[0,0],[1],[1]]:
					#print start,end
				start=end				
			for nodes in nodes_matched:
				#print nodes	
				newnetwork.node[nodes]['deg']-=newgraph.node[newmap[nodes]]['deg']
				if newnetwork.node[nodes]['deg']==0:
					newnetwork.remove_node(nodes)
				elif newnetwork.node[nodes]['deg']<0:
					raise Exception, "negative value for degeneracy after subtraction. possible error in algorithm. please try debugging"
			particle_content_list.append(highest_list)
			#index+=1
		#particle_content[index]=particle_content_list
		print "finding singlet states and higgs candidates"
		singlet_list=[]
		for nodes in newnetwork_store:
			try:
				a=newnetwork_store.node[nodes]['deg']-len(newnetwork_store.node[nodes]['vectors'])
				if a>0:
					#print newnetwork_store.node[nodes]['deg'],len(newnetwork_store.node[nodes]['vectors'])
					#print "its a"
					act_node=newnetwork_store.node[nodes]['store']
					span=newnetwork_store.node[nodes]['deg']
					indices=[]
					for states in states_dict:
						if states_dict[states][1]==act_node:
							indices.append(states-1)
					#mem=[]
					result=[]
					#print "indices are", indices
					#print "span is", span
					for i in range(span):
						newvec=form_col_vec_ind(states_dict,indices,i)
						ans=gram_schmidt_orthogonal(newvec,newnetwork_store.node[nodes]['vectors'])
						if ans!=None:
							result.append(ans)
							newnetwork_store.node[nodes]['vectors'].append(normalize(ans))
					#print "size is", len(mem)
					#print "a is", a
					
					#result=gram_schmidt_orthogonalize_many(newnetwork_store.node[nodes]['vectors'],mem[:a])
					#print result
					singlet=[]
					for elements in highest_list:
						subsinglet=[]
						for i in range(len(elements)):
							subsinglet.append(0)
						singlet.append(subsinglet)
					for i in range(len(result)):
						particle_content[index].append([singlet,[normalize(result[i])]])
						singlet_list.append(normalize(result[i]))
						#print "normal is", normalize(result[i])	
				#elif a<0:
					#print newnetwork_store.node[nodes]['deg'],len(newnetwork_store.node[nodes]['vectors'])
					#print "its not a"
			except KeyError:
				#print "its keyerror"
				act_node=newnetwork_store.node[nodes]['store']
				#print "act_node is", act_node
				#print newnetwork_store.node[nodes]['vectors']
				indices=[]
				for states in states_dict:
					if states_dict[states][1]==act_node:
						indices.append(states-1)
				result=[]
				for i in range(len(indices)):
					result.append(form_col_vec_ind(states_dict,indices,i))
				singlet=[]
				for elements in highest_list:
					subsinglet=[]
					for i in range(len(elements)):
						subsinglet.append(0)
					singlet.append(subsinglet)
				for i in range(len(result)):
					particle_content[index].append([singlet,[result[i]]])
					singlet_list.append(result[i])
					
			# Does this representation break this step? If yes, what is the state that should get the vev for this step of breaking?
		helps=0
		#print "len is", len(singlet_list)
		if index>0:
			for vec in singlet_list:
			#print "vec is", vec
			#print "done"
			#print "particle is", particle_content[index]
				for j in particle_content[index-1]:
					newlist=j[1]
					#print newlist
				#print vec,newlist
					res=gram_schmidt_orthogonal(vec,newlist)
				#print "res is", res, "done"
					if res==None:
						representation=j[0]
					#print "rep is", j[0]
						for elements in range(len(break_info_parent[index-1])):
							if break_info_parent[index-1][elements]==1:
								broken=elements
								break
						for things in representation[broken]:
							if things!=0:
								helps=1
								try:
									higgs_info[index].append(vec)
								except KeyError:
									higgs_info[index]=[vec]
								break
						break
			if helps==0:
				higgs_info[index]=[None]
		else:
			higgs_info[index]=[None]
		#for i in particle_content[index]:
			#print i[0],len(i[1])
							
						
					
						
	
							
					
						
		#			
		#			tup=generate_trial_tuples(newnetwork_store.node[nodes]['deg'],a)
		#			sol_set=[]
		#			for elements in tup:
		#				while len(sol_set)<a:
		#					matrix=form_sing_mat(elements,newnetwork_store.node[nodes]['vectors'],i,j)
		#					col_vec=form_sing_col(elements,newnetwork_store.node[nodes]['vectors'],i,j)
		#					try:
		#						B=matrix.LUsolve(col_vec)
		#						B=normalize(B)
		#						check=check_li(sol_set,B)
		#						if check=='yes':
		#							sol_set.append(B)
		#							newnetwork_store.node[nodes]['vectors'].append(B)
		#					except ValueError or ZeroDivisionError:
		#						pass
		#			if len(sol_set)<a:
		#				return Exception, "linear combination for singlet state not found"
		#			else:
		#				new_sol_set=[]
		#				for i in range(len(sol_set)):
		#					new_sol_set.append(form_sing_col_vec(sol_set[i]),indices)
		#				singlet=[]
		#				for elements in highest_list:
		#					subsinglet=[]
		#					for i in range(len(elements)):
		#						subsinglet.append(0)
		#					singlet.append(subsinglet)
		#				for i in range(sol_set):
		#					particle_content[index].append([singlet,[sol_set[i]]])
		#					
							
						 							
									
								
	return [particle_content,higgs_info,U1_gen,gen_dict]
	#return 1	

def break_rep_su2(family,dimension,highest_weight,ad):
	global total
	total=dimension
	sroots=find_sroots(family,dimension)
	output=find_weights(family,dimension,highest_weight,ad)
	#p.dump(output,g)
	#g.close()
	#output=p.load(g)	
	info=find_proots(family,dimension,sroots)
	proots=info[0]
	com_dict=info[2]
	pdim=len(info[0])
	states_dict=output[1]
	scalar_product_dict=output[2]
	highest_weight_vector=output[3]
	dimension_rep=output[4]
	mat=find_matrices(dimension,sroots,proots,com_dict,states_dict,scalar_product_dict,highest_weight_vector,dimension_rep,pdim)
	U1_gen={1:mat[0]}
	U1_vectors={1:[1,]}
	particle_content={}
	statelist=[]
	for i in range(len(states_dict)):
		statelist.append(form_col_vec_ind(states_dict,range(len(states_dict)),i))
	newlist=[[highest_weight],statelist]
	particle_content[0]=[newlist]
	particle_content[1]=[]
	higgs_info={}
	higgs_info[1]=[]
	for i in particle_content[0][0][1]:
		particle_content[1].append([[],[i]])
		higgs_info[1].append(i)
	return [particle_content,higgs_info,U1_gen]
	
		
def convert(namelist,numlist,parlist,genlist,U1list,higgs):
	l=[]
	genlist1=[]
	for things in genlist:
		genlist1.append({})
	U1list1={}
	higgs1={}
	for keys in higgs:
		higgs1[keys]=[]
	for i in parlist:
		#print i[0],len(i[1])
		l+=i[1]
	dim_rep=len(l)
	#for i in range(len(l)):
	#	for j in range(len(l)):
	#		if i!=j:
	#			prod=(l[i].T)*l[j]
	#			if prod[0]!=0:
	#				#print i,j,prod
	m=[]
	for i in range(dim_rep): 
		newvec=[0 for z in range(dim_rep)]
		for elem in range(l[i].shape[0]):
			if l[i][elem]!=0:
				newvec[elem]=l[i][elem]
		m.append(newvec)
	matrix=Matrix(m)
	#print matrix
	matrix=matrix.T
	for i in range(matrix.shape[0]):
		for j in range(matrix.shape[1]):
			matrix[i*matrix.shape[1]+j]=together(matrix[i*matrix.shape[1]+j])
	matrix1=matrix.inv()
	for i in range(matrix1.shape[0]):
		for j in range(matrix1.shape[1]):
			matrix1[i*matrix1.shape[1]+j]=together(matrix1[i*matrix1.shape[1]+j])
	#print matrix1
	for i in range(len(genlist)):
		for keys in genlist[i]:
			genlist1[i][keys]=(matrix.T)*genlist[i][keys]*matrix
	for keys in U1list:
		U1list1[keys]=(matrix.T)*U1list[keys]*matrix
	for keys in higgs:
		if higgs[keys]!=[None]:
			#for i in range(len(higgs[keys])):
			#	for j in range(len(l)):
			#		if higgs[keys][i]==l[j]:
			#			m=[0 for k in range(dim_rep)]
			#			m[j]=1
			#			higgs1[keys].append(Matrix(m))
			#			break 
				#higgs[keys][i]=matrix1*higgs[keys][i]
			for i in range(len(higgs[keys])):
				higgs1[keys].append(matrix1*higgs[keys][i])
	return [namelist,numlist,genlist1,U1list1,higgs1]
	
def convert_su2(numlist,parlist,U1list,higgs):
	l=[]
	U1list1={}
	higgs1={}
	for keys in higgs:
		higgs1[keys]=[]
	for i in parlist:
		l+=i[1]
	dim_rep=len(l)
	m=[]
	for i in range(dim_rep): 
		newvec=[0 for z in range(dim_rep)]
		for elem in range(l[i].shape[0]):
			if l[i][elem]!=0:
				newvec[elem]=l[i][elem]
		m.append(newvec)
	matrix=Matrix(m)
	matrix=matrix.T
	#print matrix.T
	for i in range(matrix.shape[0]):
		for j in range(matrix.shape[1]):
			matrix[i*matrix.shape[1]+j]=together(matrix[i*matrix.shape[1]+j])
	matrix1=matrix.inv()
	for i in range(matrix1.shape[0]):
		for j in range(matrix1.shape[1]):
			matrix1[i*matrix1.shape[1]+j]=together(matrix1[i*matrix1.shape[1]+j])
	#print matrix1
	for keys in U1list:
		U1list1[keys]=(matrix.T)*U1list[keys]*matrix
	for keys in higgs:
		if higgs[keys]!=[None]:
			for i in range(len(higgs[keys])):
				higgs1[keys].append(matrix1*higgs[keys][i])
	return [numlist,U1list1,higgs1]	


def broken_product_space(orig_dim_list,matlist,U1_gen_dict,higgsdict):
	for i in range(len(matlist)):
		for keys in matlist[i]:
			matlist[i][keys]=tens_prod(matlist[i][keys],orig_dim_list,i)
	for keys in U1_gen_dict:
		if type(keys)==tuple:
			index=keys[0]
			U1_gen_dict[keys]=tens_prod(U1_gen_dict[keys],orig_dim_list,index)
		else:
			U1_gen_dict[keys]=tens_prod(U1_gen_dict[keys],orig_dim_list,keys)
	for keys in higgsdict:
		index=keys[0]
		for j in range(len(higgsdict[keys])):
			if higgsdict[keys][j]!=None:
				higgsdict[keys][j]=tens_prod_higgs(higgsdict[keys][j],orig_dim_list,index)
	return [matlist,U1_gen_dict,higgsdict]
							
	
