from grp_rep_func import * 
from breaking import *


class LieGroup:

	def __init__(self,*args):
		if args[0]=='u' and args[1]==1:
			self.fam='u'
			self.dim=1
			self.charge=args[2]

		else:
			info=change_notation(args[0],args[1])
			self.fam=info[0]
			self.dim=info[1]
			self.broken=0

	def sroots(self):
		if self.fam=='u' and self.dim==1:
			print "Warning: The LieGroup U(1) has no simple roots"
			return []
		elif hasattr(self,'sroots_mem'):
			return self.sroots_mem
		else:
			self.sroots_mem=find_sroots(self.fam,self.dim)
			return self.sroots_mem 

	def proots(self):
		if self.fam=='u' and self.dim==1:
			print "Warning: The LieGroup U(1) has no positive roots"
			return []
		elif hasattr(self,'proots_mem'):
			return self.proots_mem
		else: 
			sroots=find_sroots(self.fam,self.dim)
			info=find_proots(self.fam,self.dim,sroots)
			#self.proots_dict=info[0]
			#l=[]
			#for keys in self.proots_dict:
			#	l.append(keys)
			self.proots_mem=info[0]
			return self.proots_mem
	
	def cmatrix(self):
		if self.fam=='u' and self.dim==1:
			print "Warning: The LieGroup U(1) has no attribute cmatrix"
			return Matrix()
		elif hasattr(self,'cartan_matrix'):
			return self.cartan_matrix
		else:
			info=find_proots(self.fam,self.dim,self.sroots())
			self.cartan_matrix=info[1]

			return Matrix(self.cartan_matrix)

	def fweights(self):
		if self.fam=='u' and self.dim==1:
			print "Warning: The LieGroup U(1) has no fundamental weights"
			return []
		elif hasattr(self,'fweights_mem'):
			return self.fweights_mem
		else:
			self.fweights_mem=find_fweights(self.dim,self.sroots())
			return self.fweights_mem

	def adim(self):
		if self.fam=='u' and self.dim==1:
			return 1
		elif hasattr(self,'adim_mem'):
			return self.adim_mem
		else:
			sroots=find_sroots(self.fam,self.dim)
			self.adim_mem=2*len(find_proots(self.fam,self.dim,sroots)[0])+len(sroots)
			return self.adim_mem
	def cdim(self):
		if self.fam=='u' and self.dim==1:
			return 1
		elif hasattr(self,'cdim_mem'):
			return self.cdim_mem
		else:
			self.cdim_mem=self.dim
			return self.cdim_mem

	def find_adjoint(self):
		if self.fam=='u' and self.dim==1:
			return None
		elif hasattr(self,'adjoint'):
			return self.adjoint
		else:
			hr=self.proots()[-1]
			res=tuple(0 for i in range(self.dim))
			sr=self.sroots()
			for i in range(len(hr)):
				res=add(res,mult(sr[i],hr[i]))
			ad=[0 for i in range(self.dim)]
			for i in range(self.dim):
				ad[i]=2*vector_product(res,sr[i])/vector_product(sr[i],sr[i])
			self.adjoint=ad
			return self.adjoint

	def convert_dl_to_root(self,dl):
		if hasattr(self,'dl_to_root_dict'):
			try:
				return self.dl_to_root_dict[dl]
			except KeyError:
				ans=self.convert_dl_to_root_pr(dl)
				self.dl_to_root_dict[dl]=ans
				return ans
		else:
			self.dl_to_root_dict={}
			ans=self.convert_dl_to_root_pr(dl)
			self.dl_to_root_dict[dl]=ans
			return ans
	
	def convert_dl_to_root_pr(self,dl):
		fw=self.fweights()
		sr=self.sroots()
		res=tuple(0 for i in range(self.dim))
		for i in range(len(dl)):
			res=add(res,mult(fw[i],dl[i]))
		root=()
		for i in range(self.dim):
			root+=(2*vector_product(res,fw[i])/vector_product(sr[i],sr[i]),)
		return root
		

	def Break(self):
		if self.fam=='u' and self.dim==1:
			print "Sorry, the Lie Group U(1) cannot be broken further"
		else:
			self.broken=0
			family=self.fam
			dimension=self.dim
			output=break_algebra([[family,dimension]])
			if output[0][0]>0:
				self.broken=1
				self.count=output[0][0]
				self.process=output[1][0]
				self.process_info=output[2][0]
				self.network_info=output[5][0]
				self.break_info_parent=output[3][0]
				self.break_info_child=output[4][0]

				# this is for automated user input to the Break function. break_history might need formatting for actual use.
				self.break_history=output[5]
			else:
				print "Warning: The algebra has not been broken" 
		#return output

class Rep:
	
	def __init__(self,liegroup,highest_weight):
		#if type(liegroup)==LieGroup:
		self.lie=liegroup
		if len(highest_weight)!=self.lie.cdim():
			if liegroup.fam=='u' and liegroup.dim==1 and len(highest_weight)==0:
				pass
			else:
				raise Exception, "Dimensionality of highest weight is incorrect"
		#else:
		#raise Exception, "first argument should be a LieGroup object"
		if self.lie.fam=='u' and self.lie.dim==1:
			if highest_weight!=[]:
				raise Exception, "The LieGroup U(1) can only have [] as its highest weight"
			else:
				self.hweight=highest_weight
				self.is_adjoint=0
		else: 
			self.hweight=highest_weight
			if self.hweight==self.lie.find_adjoint():
				self.is_adjoint=1
			else:
				self.is_adjoint=0
		
	def weights(self):
		lie=self.lie
		if lie.fam=='u' and lie.dim==1:
			print "Warning: The LieGroup U(1) does not have the attribute weights" 
			return []
		elif hasattr(self,'weights_mem'):
			return self.weights_mem
		else:
			
			hweight=self.hweight
			ad=self.is_adjoint
			sroots=self.lie.sroots()
			output=find_weights(lie.fam,lie.dim,hweight,ad)
			l=[]
			for keys in output[1]:
				l.append(output[1][keys][1])
			self.weights_mem=l
			return self.weights_mem 
	
	def dim(self):
		lie=self.lie
		if lie.fam=='u' and lie.dim==1:
			return 1 
		elif hasattr(self,'dim_mem'):
			return self.dim_mem
		else:
			
			hweight=self.hweight
			ad=self.is_adjoint
			output=find_weights(lie.fam,lie.dim,hweight,ad)
			self.dim_mem=output[4]
			return self.dim_mem
	
	def mat(self):
		lie=self.lie
		if lie.fam=='u' and lie.dim==1:
			return {0:Matrix([[lie.charge]])}
		elif hasattr(self,'mat_mem'):
			return self.mat_mem
		else:
			hweight=self.hweight
			ad=self.is_adjoint
			sroots=find_sroots(lie.fam,lie.dim)
			proots=find_proots(lie.fam,lie.dim,sroots)[0]
			out=find_weights(lie.fam,lie.dim,hweight,ad)
			info=find_proots(lie.fam,lie.dim,sroots)
			com_dict=info[2]
			pdim=len(info[0])
			states_dict=out[1]
			scalar_product_dict=out[2]
			highest_weight_vector=out[3]
			dimension_rep=out[4]
			output=find_matrices(lie.dim,sroots,proots,com_dict,states_dict,scalar_product_dict,highest_weight_vector,dimension_rep,pdim)
			self.mat_mem=output
			return self.mat_mem 

	def Break(self):
		lie=self.lie
		self.broken=0
		if lie.fam=='u' and lie.dim==1:
			print "This is a representation corresponding to an U(1) algebra. The Lie Group U(1) cannot be broken further"
		elif (lie.fam=='a' and lie.dim==1):
			if (hasattr(lie,'broken') and lie.broken==1):
				hweight=self.hweight
				ad=self.is_adjoint
				output=break_rep_su2(lie.fam,lie.dim,hweight,ad)
				self.particle_content=output[0]
				self.higgs_info=output[1]
				self.U1_gen=output[2]
				self.broken=1
			else:
				raise Exception, "Break the algebra first!"
		else:
			if (hasattr(lie,'broken') and lie.broken==1):
				count=lie.count
				process=lie.process
				process_info=lie.process_info
				network_info=lie.network_info
				break_info_parent=lie.break_info_parent
				break_info_child=lie.break_info_child
				ad=self.is_adjoint
				output=break_rep(lie.fam,lie.dim,self.hweight,count,process,process_info,break_info_parent,break_info_child,network_info,ad)
				self.particle_content=output[0]
				self.higgs_info=output[1]
				self.U1_gen=output[2]
				self.gen_list=output[3]
				self.rem_gen_list=output[4]
				self.broken=1
			else:
				raise Exception, "Break the algebra first!"	

	def to_broken_basis(self):
		liegroup=self.lie
		if hasattr(self,'broken') and self.broken==1:
			if liegroup.fam=='a' and liegroup.dim==1:
				namelist=[[]]
				numlist=[len(self.particle_content[liegroup.count][i][1]) for i in range(len(self.particle_content[liegroup.count]))]
				translist=[self.particle_content[liegroup.count][i][0] for i in range(len(self.particle_content[liegroup.count]))]
				parlist=self.particle_content[liegroup.count]
				U1list=self.U1_gen
				higgs=self.higgs_info
				output=convert_su2(numlist,parlist,U1list,higgs)
				self.broken_rep_dim=output[0]
				self.broken_U1_gen=output[1]
				self.higgs=output[2]
				self.broken_group_list=namelist
				self.trans=translist
				self.changed=1
			else:
				namelist=liegroup.process_info[liegroup.count]
				numlist=[len(self.particle_content[liegroup.count][i][1]) for i in range(len(self.particle_content[liegroup.count]))]
				translist=[self.particle_content[liegroup.count][i][0] for i in range(len(self.particle_content[liegroup.count]))]
				parlist=self.particle_content[liegroup.count]
				genlist=self.gen_list[liegroup.count]
				rem_genlist=self.rem_gen_list[liegroup.count]
				U1list=self.U1_gen
				higgs=self.higgs_info
				output=convert(namelist,numlist,parlist,genlist,rem_genlist,U1list,higgs)
				self.broken_group_list=output[0]
				self.broken_rep_dim=output[1]
				self.broken_gen_list=output[2]
				self.broken_rem_gen_list=output[3]
				self.broken_U1_gen=output[4]
				self.higgs=output[5]
				self.trans=translist
				self.changed=1
		else:
			raise Exception, "The representation has not been broken. Please break it first"
					
					 
			


class ProdGroup:
	
	def __init__(self,list_of_groups,list_of_U1grp):
		for grp in list_of_groups:
			if grp.fam=='u'and grp.dim==1:
				return Exception, "Enter U1 groups separately"
		self.list=list_of_groups
		self.U1list=list_of_U1grp

	def Break(self):
		temp=[]
		for each in self.list:
			temp.append([each.fam,each.dim])
		output=break_algebra(temp)
		for index in range(len(self.list)):
			if output[0][index]>0:
				self.list[index].broken=1
				self.list[index].count=output[0][index]
				self.list[index].process=output[1][index]
				self.list[index].process_info=output[2][index]
				self.list[index].network_info=output[5][index]
				self.list[index].break_info_parent=output[3][index]
				self.list[index].break_info_child=output[4][index]

				# this is for automated user input to the Break function. break_history might need formatting for actual use.
				self.list[index].break_history=output[5][index]
					
	
	def sroots(self,num):
		try:
			return self.list[num].sroots()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"

	def proots(self,num):
		try:
			return self.list[num].proots()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"

	def fweights(self,num):
		try:
			return self.list[num].fweights()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"

	def cmatrix(self,num):
		try:
			return self.list[num].cmatrix()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"

	def adim(self,num):
		try:
			return self.list[num].adim()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"

	def cdim(self,num):
		try:
			return self.list[num].cdim()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"


class ProdGroupRep:
	
	def __init__(self,prodgroup,list_of_highest_weights):
		lielist=prodgroup.list
		if len(lielist)!=len(list_of_highest_weights):
			raise Exception, "Dimension of list of highest weights does not match that of list of lie groups"
		l=[]
		for  i in range(len(lielist)):
			l.append(Rep(lielist[i],list_of_highest_weights[i]))
		self.list=l

	
	def weights(self,num):
		try:
			return self.list[num].weights()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"
		
	def dim(self,num):
		try:
			return self.list[num].dim()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"	
				
	def mat(self,num):
		try:
			return self.list[num].mat()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"	
		
	def Break(self,num):
		try:
			self.list[num].Break()
		except IndexError:
			raise Exception, "Sorry, cannot break the group since there are no LieGroups corresponding to this index"

	def Breakall(self):
		for each in self.list:
			if hasattr(each.lie,'broken') and each.lie.broken==1:
				each.Break()		

	
	#def summarize(self)	

	# This function might not be required anymore. Anyway, it has to be modified because of the introduction of the separate U1list in ProdGroups	

	#def to_broken_prod_basis(self):
	#	replist=self.list
	#	l=[]
	#	for i in range(len(replist)):
	#		rep=replist[i]
	#		liegroup=replist[i].lie	
	#		if hasattr(liegroup,'broken') and liegroup.broken==1:
	#			if hasattr(rep,'broken') and rep.broken==1:
	#				if liegroup.fam=='a' and liegroup.dim==1:
	#					if hasattr(rep,'changed') and rep.changed==1:
	#						l.append(2)
	#					else:
	#						rep.to_broken_basis()
	#						l.append(2)
	#				else:
	#					if hasattr(rep,'changed') and rep.changed==1:
	#						l.append(1)
	#					else:
	#						rep.to_broken_basis()
	#						l.append(1)
	#			else:
	#				while True:
	#					inp=raw_input("This rep is not broken. Do you want to break it? yes/no")
	#					if inp=='yes':
	#						rep.Break()
	#						rep.to_broken_basis()
	#						if liegroup.fam=='a' and liegroup.dim==1:
	#							l.append(2)
	#							break
	#						else:
	#							l.append(1)
	#							break
	#					elif inp=='no':
	#						l.append(0)
	#						break
	#					else:
	#						print "Invalid input. Try again"
	#		elif liegroup.fam=='u' and liegroup.dim==1:
	#			l.append(3)
	#		else:
	#			l.append(0)
	#	grouplist=[]
	#	dimlist=[]
	#	orig_dim_list=[]
	#	matlist=[]
	#	#parlist=[]
	#	translist=[]
	#	U1_gen_dict={}
	#	higgsdict={}
	#	for i in range(len(l)):
	#		if l[i]==0:
	#			dim=replist[i].dim()
	#			grouplist.append([replist[i].lie.fam,replist[i].lie.dim])
	#			dimlist.append([dim])
	#			orig_dim_list.append(dim)
	#			matlist.append(replist[i].mat())
	#			translist.append([replist[i].hweight])
	#		elif l[i]==1:
	#			grouplist+=replist[i].broken_group_list
	#			dimlist.append(replist[i].broken_rep_dim)
	#			orig_dim_list.append(replist[i].dim())
	#			matlist+=replist[i].broken_gen_list
	#			for keys in replist[i].broken_U1_gen:
	#				tup=(i,keys)
	#				U1_gen_dict[tup]=replist[i].broken_U1_gen[keys]
	#			for keys in replist[i].higgs:
	#				tup=(i,keys)
	#				higgsdict[tup]=replist[i].higgs[keys]
	#			translist.append(replist[i].trans)
	#		elif l[i]==2:
	#			grouplist+=replist[i].broken_group_list
	#			dimlist.append(replist[i].broken_rep_dim)
	#			orig_dim_list.append(replist[i].dim())
	#			matlist.append({})
	#			for keys in replist[i].broken_U1_gen:
	#				tup=(i,keys)
	#				U1_gen_dict[tup]=replist[i].broken_U1_gen[keys]
	#			for keys in replist[i].higgs:
	#				tup=(i,keys)
	#				higgsdict[tup]=replist[i].higgs[keys]
	#			translist.append(replist[i].trans)
	#		elif l[i]==3:
	#			grouplist.append([[]])
	#			dimlist.append([1])
	#			orig_dim_list.append(replist[i].dim())
	#			matlist.append({})
	#			U1_gen_dict[i]=replist[i].mat()[0]
	#			translist.append([])
	#	#print higgsdict
	#	output=broken_product_space(orig_dim_list,matlist,U1_gen_dict,higgsdict)
	#	self.broken_group_list=grouplist
	#	self.broken_rep_dim=dimlist
	#	self.trans=translist
	#	self.broken_U1_gen=output[1]
	#	self.broken_gen_list=output[0]
	#	self.higgs=output[2]
	#	self.changed=1
				
					
							
		
		
	
	
