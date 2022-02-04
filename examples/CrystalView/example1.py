from XtraCrysPy import XtraCrysPy as XCP

if '__main__' == __name__:

  # Define params
  a = 10.599478

  # Unit origin
  origin = [0,0,0]

  # Cell parameters (Normalized)
  fcc_vecs = [[a/2,a/2,0], [0,a/2,a/2], [a/2,0,a/2]]

  # Atomic basis
  atoms = [[0,0,0], [a/4,a/4,a/4]]

  # Label for each atom and colors
  labels = ['Ga', 'As']

  # Bond distance between Ga and As, in Anstrom
  bonds = {'Ga_As':8}

  # Species information, including radius and color
  species = {'Ga':{'color':(1,0,0,0),'radius':.8}, 'As':{'color':(0,0,1,1)}}

  # Construct model with XtraCrysPy
  xcpy = XCP.XtraCrysPy(lattice=fcc_vecs, basis=atoms, species=species, basis_labels=labels, origin=origin, bonds=bonds)

  # Write an XML file for blender
  xcpy.write_blender_xml(fname='GaAs_blender.xml', nx=3, ny=3, nz=3)

  xcpy.start_cryspy_view(title='GaAs', nx=3, ny=3, nz=3)

