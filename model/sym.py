from sympy import *

def CreateSymbols(group,name,*args):
	if group=='f':
		#the order of arguments is particle_type, particle_index, family_index, gauge_group, representation, matter_index, barred and chirality
		a=Fermion(name,commutative=False)	

		# 1 denotes either modification or definition and is followed by an argument, 0 denotes display for these attributes. 
		a.part_type(1,'f')
		a.part_ind(1,args[0])
		a.fam_ind(1,args[1])
		a.ggrp(1,args[2])
		a.rep(1,args[3])
		a.mat_ind(1,args[4])
		a.barred(1,args[5])
		a.chirality(1,args[6])
		return a

	if group=='gs':
		# the order of arguments is type,gauge_group,representation,matter_rep,dynkin labels and origin 

		a=GaugeS(name,commutative=False)
		a.gtype(1,'gs')
		a.ggrp(1,args[0])
		a.rep(1,args[1])
		a.mat_rep(1,args[2])
		a.dl(1,args[3])
		a.ori(1,args[4])
		return a

	if group=='gl':
		# the order of arguments is type,gauge_group,representation,matter_rep,dynkin labels,origin and lorentz

		a=GaugeL(name,commutative=False)
		a.gtype(1,'gl')
		a.ggrp(1,args[0])
		a.rep(1,args[1])
		#a.mat_rep(1,args[3])
		a.dl(1,args[2])
		a.ori(1,args[3])
		a.lorentz(1,args[4])
		a.pos(1,args[5])	
		return a


class Fermion(Symbol):
	def part_type(self,*args):
		if args[0]==0:
			if hasattr(self,'part_type'):
				return self.part_type
			else:
				return Exception, "the fermionic field %s does not have the attribute part_type" % self.name
		elif args[1]==1:
			self.part_type=args[1]
		else:
			raise Exception, "the constructor part_type takes only 1 (write) or 0 (display) as its first argument"

	def part_ind(self,*args):
		if args[0]==0:
			if hasattr(self,'part_ind'):
				return self.part_ind
			else:
				return Exception, "the fermionic field %s does not have the attribute part_ind" % self.name
		elif args[1]==1:
			self.part_ind=args[1]
		else:
			raise Exception, "the constructor part_ind takes only 1 (write) or 0 (display) as its first argument"
	def fam_ind(self,*args):
		if args[0]==0:
			if hasattr(self,'fam_ind'):
				return self.fam_ind
			else:
				return Exception, "the fermionic field %s does not have the attribute fam_ind" % self.name
		elif args[1]==1:
			self.fam_ind=args[1]
		else:
			raise Exception, "the constructor fam_ind takes only 1 (write) or 0 (display) as its first argument"	
	def ggrp(self,*args):
		if args[0]==0:
			if hasattr(self,'ggrp'):
				return self.ggrp
			else:
				return Exception, "the fermionic field %s does not have the attribute ggrp" % self.name
		elif args[1]==1:
			self.ggrp=args[1]
		else:
			raise Exception, "the constructor ggrp takes only 1 (write) or 0 (display) as its first argument"
	def rep(self,*args):
		if args[0]==0:
			if hasattr(self,'rep'):
				return self.rep
			else:
				return Exception, "the fermionic field %s does not have the attribute rep" % self.name
		elif args[1]==1:
			self.rep=args[1]
		else:
			raise Exception, "the constructor rep takes only 1 (write) or 0 (display) as its first argument"
	def mat_ind(self,*args):
		if args[0]==0:
			if hasattr(self,'mat_ind'):
				return self.mat_ind
			else:
				return Exception, "the fermionic field %s does not have the attribute mat_ind" % self.name
		elif args[1]==1:
			self.mat_ind=args[1]
		else:
			raise Exception, "the constructor mat_ind takes only 1 (write) or 0 (display) as its first argument"
	def barred(self,*args):
		if args[0]==0:
			if hasattr(self,'barred'):
				return self.barred
			else:
				return Exception, "the fermionic field %s does not have the attribute barred" % self.name
		elif args[1]==1:
			self.barred=args[1]
		else:
			raise Exception, "the constructor barred takes only 1 (write) or 0 (display) as its first argument"
	def chirality(self,*args):
		if args[0]==0:
			if hasattr(self,'chirality'):
				return self.chirality
			else:
				return Exception, "the fermionic field %s does not have the attribute chirality" % self.name
		elif args[1]==1:
			self.chirality=args[1]
		else:
			raise Exception, "the constructor chirality takes only 1 (write) or 0 (display) as its first argument"


class GaugeS(Symbol):
	def gtype(self,*args):
		if args[0]==0:
			if hasattr(self,'gtype'):
				return self.gtype
			else:
				return Exception, "the GaugeS field %s does not have the attribute gtype" % self.name
		elif args[0]==1:
			self.gtype=args[1]
		else:
			raise Exception, "the constructor gtype takes only 1 (write) or 0 (display) as its first argument"
	def ggrp(self,*args):
		if args[0]==0:
			if hasattr(self,'ggrp'):
				return self.ggrp
			else:
				return Exception, "the GaugeS field %s does not have the attribute ggrp" % self.name
		elif args[0]==1:
			self.ggrp=args[1]
		else:
			raise Exception, "the constructor ggrp takes only 1 (write) or 0 (display) as its first argument"
	def rep(self,*args):
		if args[0]==0:
			if hasattr(self,'rep'):
				return self.rep
			else:
				return Exception, "the GaugeS field %s does not have the attribute rep" % self.name
		elif args[0]==1:
			self.rep=args[1]
		else:
			raise Exception, "the constructor rep takes only 1 (write) or 0 (display) as its first argument"
	def mat_rep(self,*args):
		if args[0]==0:
			if hasattr(self,'mat_rep'):
				return self.mat_rep
			else:
				return Exception, "the GaugeS field %s does not have the attribute mat_rep" % self.name
		elif args[0]==1:
			self.mat_rep=args[1]
		else:
			raise Exception, "the constructor mat_rep takes only 1 (write) or 0 (display) as its first argument"
	def dl(self,*args):
		if args[0]==0:
			if hasattr(self,'dl'):
				return self.dl
			else:
				return Exception, "the GaugeS field %s does not have the attribute dl" % self.name
		elif args[0]==1:
			self.dl=args[1]
		else:
			raise Exception, "the constructor dl takes only 1 (write) or 0 (display) as its first argument"
	def ori(self,*args):
		if args[0]==0:
			if hasattr(self,'ori'):
				return self.ori
			else:
				return Exception, "the GaugeS field %s does not have the attribute ori" % self.name
		elif args[0]==1:
			self.ori=args[1]
		else:
			raise Exception, "the constructor ori takes only 1 (write) or 0 (display) as its first argument"
	
class GaugeL(Symbol):
	def gtype(self,*args):
		if args[0]==0:
			if hasattr(self,'gtype'):
				return self.gtype
			else:
				return Exception, "the GaugeL field %s does not have the attribute gtype" % self.name
		elif args[0]==1:
			self.gtype=args[1]
		else:
			raise Exception, "the constructor gtype takes only 1 (write) or 0 (display) as its first argument"
	def ggrp(self,*args):
		if args[0]==0:
			if hasattr(self,'ggrp'):
				return self.ggrp
			else:
				return Exception, "the GaugeL field %s does not have the attribute ggrp" % self.name
		elif args[0]==1:
			self.ggrp=args[1]
		else:
			raise Exception, "the constructor ggrp takes only 1 (write) or 0 (display) as its first argument"
	def rep(self,*args):
		if args[0]==0:
			if hasattr(self,'rep'):
				return self.rep
			else:
				return Exception, "the GaugeL field %s does not have the attribute rep" % self.name
		elif args[0]==1:
			self.rep=args[1]
		else:
			raise Exception, "the constructor rep takes only 1 (write) or 0 (display) as its first argument"
	def pos(self,*args):
		if args[0]==0:
			if hasattr(self,'mat_rep'):
				return self.mat_rep
			else:
				return Exception, "the GaugeL field %s does not have the attribute mat_rep" % self.name
		elif args[0]==1:
			self.pos=args[1]
		else:
			raise Exception, "the constructor mat_rep takes only 1 (write) or 0 (display) as its first argument"
	def dl(self,*args):
		if args[0]==0:
			if hasattr(self,'dl'):
				return self.dl
			else:
				return Exception, "the GaugeL field %s does not have the attribute dl" % self.name
		elif args[0]==1:
			self.dl=args[1]
		else:
			raise Exception, "the constructor dl takes only 1 (write) or 0 (display) as its first argument"
	def ori(self,*args):
		if args[0]==0:
			if hasattr(self,'ori'):
				return self.ori
			else:
				return Exception, "the GaugeL field %s does not have the attribute ori" % self.name
		elif args[0]==1:
			self.ori=args[1]
		else:
			raise Exception, "the constructor ori takes only 1 (write) or 0 (display) as its first argument"
	def lorentz(self,*args):
		if args[0]==0:
			if hasattr(self,'lorentz'):
				return self.lorentz
			else:
				return Exception, "the GaugeL field %s does not have the attribute lorentz" % self.name
		elif args[0]==1:
			self.lorentz=args[1]
		else:
			raise Exception, "the constructor lorentz takes only 1 (write) or 0 (display) as its first argument"

def CreateTrace(expr):
	#print str(expr)
	name="TR(%s)" % expr
	#print name
	res=Trace(name)
	#print res
	res.set_expr(expr)
	return res



class Trace(Symbol):
	def set_expr(self,expr):
		self.expr=expr


def CreateDS(expr):
	name='DS(%s)' % expr
	res=DS(name)
	res.set_expr(expr)
	return res



class DS(Symbol):
	def set_expr(self,expr):
		self.expr=expr

def CreateDL(expr,lorentz,pos):
	name='DL_%s_%s(%s)' % (lorentz,pos,expr)
	res=DL(name)
	res.set_expr(expr)
	res.set_lorentz(lorentz,pos)
	return res



class DL(Symbol):
	def set_expr(self,expr):
		self.expr=expr
	def set_lorentz(self,lorentz,pos):
		self.lorentz=lorentz
		self.pos=pos	







	
