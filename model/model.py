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
								print self.ggrp.list[index]
								self.broken_glag[0][index]=self.pr_break_glag(index)
						else:
							self.ggrp.list[index].Break()
							if self.ggrp.list[index].broken==1:
								self.list_broken[index]=1
								self.broken_glag[0][index]=self.pr_break_glag(index)
									
									



								
				
				
		else:
			return Exception, "define the gauge group symmetry for this model first and then try breaking the symmetry" 
      





	def pr_break_glag(self,index):
		if self.ggrp.list[index].broken==0:
			return 'intact'
		else:
			return 'broken'
			
			
					
			
			
			
			
		
			
		
