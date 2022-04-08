from XtraCrysPy import XCP_Atoms as XCP
from sys import argv

if '__main__' == __name__:

  argc = len(argv)
  if argc < 2 or argv[1] == '-h':
    print('Usage:')
    print('  python main.py <inputfile>')
    print('  python main.py <inputfile> <bond_length>')
    exit()

  b_len = 6.1
  if argc > 2:
    b_len = float(argv[2])

  fname = argv[1]
  a_info = {'bonds':b_len}
  xcp = XCP.XCP_Atoms(params=a_info, model=fname, sel_type='Chain')
  xcp.start_crystal_view()

