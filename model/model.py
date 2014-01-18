# importing the LieAlgebra module

import sys
sys.path.insert(0, '/home/dibya/adding_final/liealgebra')
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
		self.matter=[]
		
	# public function to change the name of the model at any stage

	def change_name(self,name):
		self.name=name

	# public function to change SUSY

	def add_susy(self,num):
		if self.glag_added==0:
			if num in [0,1,2,4]:
				self.pr_add_susy(num)
			else:
				return Exception, "This package has support for N=0,1,2 and 4 SUSY"
		else:
			print "Cannot modify susy at the current stage. Please start with a new model"

	# private function to change SUSY

	def pr_add_susy(self,num):
		self.susy=num
		self.susy_added=1

	# public function to change Gauge Group

	def add_ggrp(self, *args):
		if self.susy_added==1:
			if self.flag_added==0:
				if len(args)>1 and len(args[1][0])==len(args[0].list) and len(args[1][1])==len(args[0].U1list):
					self.coupling=args[1]
				else:
					coupling=[[],[]]
					for num in range(len(args[0].list)):
						coupling[0].append(Symbol('g%s' % num,commutative=False))
					for num in range(len(args[0].list),len(args[0].list)+len(args[0].U1list),1):
						coupling[1].append(Symbol('g%s' % num,commutative=False))
					self.coupling=coupling
				self.pr_add_ggrp(args[0],self.susy)
			else:
				"Cannot change gauge group anymore"
		else:
			print "Please define SUSY before entering gauge group"

	# private function to change Gauge Group

	def pr_add_ggrp(self,prdgrp,susy):
		self.ggrp=prdgrp
		self.pr_add_glag(self.ggrp,susy)

	def pr_add_glag(self,ggrp,susy):
		expr_na_gauge=[]
		expr_a_gauge=[]
		namelist=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
		glist=ggrp.list
		U1list=ggrp.U1list
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
					newsym=CreateSymbols('gl',letter+'%s_%s_0' % (lorentz_index,pos), ggrp, rep,statedict,origin,lorentz_index,pos)
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
					newsym=CreateSymbols('gl',letter+'%s_%s_0' % (lorentz_index,pos), ggrp, rep,None,origin,lorentz_index,pos)
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
		 




			
			
					
			
			
			
			
		
			
		
