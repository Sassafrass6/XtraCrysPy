from XtraCrysPy import XtraCrysPy as XCP

if '__main__' == __name__:

  species = {'Ba':{'color':(.1,1,.1), 'radius':1.4}, 'Cu':{'color':(.8,.8,0)}, 'As':{'color':(0,0,.7)}}
  bonds = {'Cu_As':5}
  crystal = XCP.XtraCrysPy(inputfile='BaCuAs.scf.in', species=species, bonds=bonds)
  crystal.start_cryspy_view()
