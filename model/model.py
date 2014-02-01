# importing the LieAlgebra module

import sys
sys.path.insert(0, '/home/dibya/K2P2/liealgebra')
from user import *
from sym import *
from sympy import *


class Model:

	# initialize model with optional name argument

	def __init__(self,*args):
		if len(args)==0:
			pass
		elif len(args)==1:
			self.name=args[0]
		else:
			return Exception, "Model object constructor can take at most 1 argument" 

		# flags

		self.susy_added=0
		self.glag_added=0
		self.flag_added=0
		self.broken=0
		self.matter=[]
		
	# public function to change the name of the model at any stage

	def change_name(self,name):
		self.name=name

	# public function to change SUSY

	def add_susy(self,num):
		if self.glag_added==0:
			if num==0:
				self.pr_add_susy(num)
			else:
				return Exception, "This package has support for N=0 currently"
		else:
			print "Cannot modify susy at the current stage. Please start with a new model"

	# private function to change SUSY

	def pr_add_susy(self,num):
		self.susy=num
		self.susy_added=1

	# public function to change Gauge Group

	def add_ggrp(self, *args):
		if self.susy_added==1:
			if self.flag_added==0 or self.broken==1:
				if len(args)>1 and len(args[1][0])==len(args[0].list) and len(args[1][1])==len(args[0].U1list):
					self.coupling=args[1]
				else:
					coupling=[[],[]]
					for num in range(len(args[0].list)):
						coupling[0].append(Symbol('g%s' % num,commutative=False))
					for num in range(len(args[0].list),len(args[0].list)+len(args[0].U1list),1):
						coupling[1].append(Symbol('g%s' % num,commutative=False))
					self.coupling=coupling
				self.pr_add_ggrp(args[0])
			else:
				"Cannot change gauge group anymore"
		else:
			print "Please define SUSY before entering gauge group"

	# private function to change Gauge Group

	def pr_add_ggrp(self,prdgrp):
		self.ggrp=prdgrp
		self.higgs_list=[{} for i in range(len(self.ggrp.list))]
		self.last_rep_limit=[None for range(len(self.ggrp.list))]
		self.pr_add_glag()

	def pr_add_glag(self):
		expr_na_gauge=[]
		expr_a_gauge=[]
		namelist=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
		glist=self.ggrp.list
		U1list=self.ggrp.U1list
		rep_init=[[],[]]
		for each in glist:
			rep_init[0].append([0,])
		for each in U1list:
			rep_init[1].append(0)
		for index in range(len(glist)):
			letter=namelist[0]
			namelist=namelist[1:]
			gtype='gl'	
			thisgrp=glist[index]
			adjoint=thisgrp.find_adjoint()
			rep=[rep_init[0][:],rep_init[1]]
			rep[0][index]=adjoint
			newrep=Rep(thisgrp,adjoint)
			matrices=newrep.mat()
			#for keys in matrices:
			#	print type(keys[0])
			#	break
			statedict={}
			cartan_counter=0
			state_counter=0
			for states in newrep.weights():
				if states!=tuple(0 for i in range(thisgrp.dim)):
					newtuple=thisgrp.convert_dl_to_root(states)
					thistuple=()
					for i in newtuple:
						thistuple+=(int(i),)
					#print type(newtuple[0])
					statedict[state_counter]=[states,matrices[thistuple]]
					state_counter+=1
				else:
					statedict[state_counter]=[states,matrices[cartan_counter]]
					cartan_counter+=1
					state_counter+=1
			origin='START'
			symbols=[]
			for lorentz_index in ['mu','nu']:
				for pos in ['up','down']:
					newsym=CreateSymbols('gl',letter+'%s_%s_0' % (lorentz_index,pos), self.ggrp, rep,statedict,origin,lorentz_index,pos)
					symbols.append(newsym)
			argument=(CreateDL(symbols[2],'mu','up')-CreateDL(symbols[0],'nu','up')+I*self.coupling[0][index]*(symbols[0]*symbols[2]-symbols[2]*symbols[0]))*(CreateDL(symbols[3],'mu','down')-CreateDL(symbols[1],'nu','down')+I*self.coupling[0][index]*(symbols[1]*symbols[3]-symbols[3]*symbols[1]))
			expr_na_gauge.append(CreateTrace(argument))	

		for index in range(len(U1list)):
			letter=namelist[0]
			namelist=namelist[1:]
			gtype='gl'
			rep=[rep_init[0][:],rep_init[1]]
			rep[1][index]=1	
			origin='START'
			symbols=[]
			for lorentz_index in ['mu','nu']:
				for pos in ['up','down']:
					newsym=CreateSymbols('gl',letter+'%s_%s_0' % (lorentz_index,pos), self.ggrp, rep,None,origin,lorentz_index,pos)
					symbols.append(newsym)
			argument=(CreateDL(symbols[2],'mu','up')-CreateDL(symbols[0],'nu','up')+I*self.coupling[1][index]*(symbols[0]*symbols[2]-symbols[2]*symbols[0]))*(CreateDL(symbols[3],'mu','down')-CreateDL(symbols[1],'nu','down')+I*self.coupling[1][index]*(symbols[1]*symbols[3]-symbols[3]*symbols[1]))
			expr_a_gauge.append(argument)
		self.glag_expr=[expr_na_gauge,expr_a_gauge]
		newexpr=0
		for i in expr_na_gauge+expr_a_gauge:
			newexpr+=i
		self.glag=newexpr
		#print self.glag
		self.glag_added=1


	def add_matter(self,*args):
		if model.glag_added==1:
				if flag_added==1:
					self.pr_add_matter(self,1,*args)
				else:	
					self.pr_add_matter(self,0,*args)
		else:
			return Exception, "Cannot add matter if gauge group is not defined"
	
	def pr_add_matter(self,bin,*args):
		if bin==0:
			self.flag=0
			for each in args:
				self.pr_add_flag(each)
				self.matter.append(each)



# I am assuming that Ritesh Sir will add the codes for adding matter to model.py. Once this is done, flag and slag will be defined and flags for functions that change gauge groups and susy will be set down. At this point the routine for breaking the algebra becomes available. This is a function which takes the index of the lie group as input. Once this function is called, an interactive breaking program will take over (this has to be automated later) asking the user about how the algebra is to be broken. Once this input is fully provided new attributes become available. These new attributes are:

#  1. particle content: tells us how the representation is broken into smaller representations and what states of the original representation correspond to smaller representations
#  2. higgs_info: tells us if there are singlets left after the breaking (Q: does this also check U1 charge? A: no it does not, it just returns singlets that break certain steps)
#  3. U1_gen: gives the U1 generator (in the new basis) formed as a result of deletion of a circle in the dynkin diagram 
#  4. gen_list: this is the list of generators of the algebra in the new broken basis(?) this currently gives only those generators that are a part of the subalgebras. Need to add codes in /liealgebra to also give matrices corresponding to other generators in the broken basis

# Based on these new attributes, we need to rewrite the lagrangian : glag, flag and slag.

# Steps

#	1. Change the relevant file in /liealgebra to give all generators in the changed basis
#	2. Add the routine for interactive breaking in model
#	3. Add private functions to change glag, flag and slag



	def Break(self, *args):
		if self.glag_added==1:

			# Breaking can proceed in the absence of matter. 
			#	1. Need to define a function later on which will be able to add matter and break the reps and write flag and slag according to the most recent breaking scheme for that model. 
			#	2. Also should add a function later on which can provide a new copy of model objects which can be independently manipulated (branches)

			# CONTENT OF THIS FUNCTION


			# Step 1: Let the user break the non abelian algebras as he pleases

			# Way 1: When user mentions which algebras he wants broken -> he calls model.Break() which arguments

			# I need to comment on this part

			self.list_broken=[0 for i in range(len(self.ggrp.list))]
				
			if self.broken==0:
				self.broken_glag=[[],self.glag_expr[1]]
				for index in range(len(self.ggrp.list)):
					if len(args)==0 or (len(args)>0 and index in args):
						self.ggrp.list[index].Break()
						if self.ggrp.list[index].broken==1:
							self.list_broken[index]=1
							self.broken_glag[0].append(self.pr_break_glag(index))
						else:
							self.broken_glag[0].append(self.pr_break_glag(index))
					else:
						self.broken_glag[0].append(self.pr_break_glag(index))
				self.broken=1
						
			else:
				for index in range(len(self.ggrp.list)):
					if len(args)==0 or (len(args)>0 and index in args):
						if self.ggrp.list[index].broken==1:
							prev_breaking=self.ggrp.list[index].break_history
							self.ggrp.list[index].Break()
							if self.ggrp.list[index].broken==0 or self.ggrp.list[index].break_history!=prev_breaking:
								self.list_broken[index]=1
								self.broken_glag[0][index]=self.pr_break_glag(index)
								self.higgs_list[index]={}
						else:
							self.ggrp.list[index].Break()
							if self.ggrp.list[index].broken==1:
								self.list_broken[index]=1
								self.broken_glag[0][index]=self.pr_break_glag(index)
								self.higgs_list[index]={}
									
									



								
				
				
		else:
			return Exception, "define the gauge group symmetry for this model first and then try breaking the symmetry" 



	def higgs_info(self, fork, rep_limit):
		# Should we get higgs_info on ProdRep s too?

		if self.broken==1
			# setting up rep_limit_dict
			if fork=='all' and type(rep_limit)==tuple:
				if len(rep_limit)=2 and type(rep_limit[0])==int and type(rep_limit)==int and rep_limit[0]>0 and rep_limit[1]>rep_limit[0]:
					rep_limit_dict={}
					for index in range(len(self.ggrp.list)):
						rep_limit_dict[index]=rep_limit
				else:
					return Exception, "invalid rep_limit"
			elif (fork=='all' or fork=='select') and type(rep_limit)==dict:
				for keys in rep_limit:
					if rep_limit[keys] in range(len(self.ggrp.list)):
						if len(rep_limit[keys])=2 and type(rep_limit[keys][0])==int and type(rep_limit[keys])==int and rep_limit[keys][0]>0 and rep_limit[keys][1]>rep_limit[keys][0]:
							pass
						else:
							return Exception, "invalid rep_limit"
					else:
						return Exception, "group index out of range"
				rep_limit_dict=rep_limit
			else:
				return Exception, "invalid input" 
				
			# checking for all groups
			for num in range(len(self.ggrp.list)):
				if num in rep_limit_dict.keys:
					limit=rep_limit_dict[num]
					# check if the group is broken
					if self.grp.list[num].broken==0:
						print "group number", num, "is not broken. cannot generate information on higgs for symmetry breaking of this group"
					# otherwise
					else:  
						# if the group was broken recently
						if self.list_broken[num]==1:
							# initializing the list. note that this should be defined at thetime of defining Model.ggrp	
							self.higgs_list[num]={}
							# getting the list out of a function 
							# the function is not yet present
							dict_of_reps=generate_reps(self.ggrp.list[num], range(limit))
							# breaking the reps and adding them to the list
							for keys in dict_of reps:
								for reps in dict_of_reps[keys]:
									rep_object=Rep(self.ggrp.list[num],reps)
									rep_object.Break()
									try:
										self.higgs_list[num][keys].append(rep_object)
									except:
										self.higgs_list[num][keys]=[rep_object,]
							# changing self.last_rep_limit to its current value
							self.last_rep_limit=limit
						# if the grp was not broken recently but rep limit has changed
						# note that self.last_rep_limit should be initialized at thetime of defining Model.ggrp
						elif self.list_broken[num]==0 and limit!=self.last_rep_limit[num]:
							# finding the new keys that need to be added
							new_keys=()
							for i in range(limit):
								if i not in self.higgs_list[num].keys:
									new_keys+=(i,)
							# removing keys that are not in the present range
							for keys in self.higgs_list[num].keys:
								if keys not in range(limit):
									del self.higgs_list[num][keys]
							# finding reps corresponding to new keys
							dict_of_reps=generate_reps(self.ggrp.list[num], range(limit))
							# breaking the reps and adding them to the list
							for keys in dict_of reps:
								for reps in dict_of_reps[keys]:
									rep_object=Rep(self.ggrp.list[num],reps)
									rep_object.Break()
									try:
										self.higgs_list[num][keys].append(rep_object)
									except:
										self.higgs_list[num][keys]=[rep_object,]
							# changing self.last_rep_limit to its current value
	 						self.last_rep_limit=limit
									
										
									 

							 
						
					
		
		
			





	def pr_break_glag(self,index):
		if self.ggrp.list[index].broken==0:
			return 'intact'
		else:
			return 'broken'
			
						
			
			
		
#---------------------------------SECONDARY FUNCTIONS--------------------------#


def generate_reps(grp,limit):
	cdim=grp.cdim()
	rep_dict={}
	for num in limit:
		rep_dict[num]=find_reps(cdim,num)

	return rep_dict

def partition(number):
	answer = set()
     	answer.add((number, ))
     	for x in range(1, number):
     	    for y in partition(number - x):
     	        answer.add(tuple(sorted((x, ) + y)))
     	return answer

def find_pos(pos_tuple, num):
	if num==0:
		answer=set()
		answer.add(tuple())
		return answer
	else:
		answer=set()
		for i in range(len(pos_tuple)):
			l=set()
			for j in find_pos(pos_tuple[:i]+pos_tuple[i+1:], num-1):
				if pos_tuple[i] not in j:
					l.add((pos_tuple[i],)+j)
			answer=answer.union(l)
		return answer
				
					
			
def find_reps(cdim,num):
	answer=[]
	for part in partition(num):
		length=len(part)
		if len(part)<=cdim:
			for pos in find_pos(tuple(j for j in range(cdim)), length):
				l=[0 for i in range(cdim)]
				ind=0
				for elements in pos:
					l[elements]=part[ind]
					ind+=1
				if l not in answer:
					answer.append(l)
	return answer
					
				  
	
	
	 
					
		
