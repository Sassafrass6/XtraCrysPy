from XtraCrysPy import XCP_Atoms

if '__main__' == __name__:

  # Label for each atom and colors
  labels = ['O']*2 + ['N']*4 + ['C']*8 + ['H']*10

  # Bond distance in Angstrom
  bonds = {'C_C':1.5, 'C_N':1.5, 'C_O':1.5, 'C_H':1.2}

  # Visual information: species, lattice, atoms, bonds, units
  a_info = {'species':labels}#, 'units':'angstrom'}

  xcp = XCP_Atoms.XCP_Atoms(model='data_files/lactic_acid.cif', params=a_info)
  xcp.start_crystal_view()

