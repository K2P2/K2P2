from user import *
grp=LieGroup('so',10)
hr=grp.proots()[-1]
hdl=grp.find_adjoint()
new=grp.convert_dl_to_root(tuple(hdl))
print new,hr
