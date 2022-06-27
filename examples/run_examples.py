'''
 Run all examples, or run examples specified by any number of input arguments.
'''

usage = '''
print('Usage:\n  python run_examples.py')
'''

if __name__ == '__main__':
  from os import system
  from sys import argv

  examples = []
  argc = len(argv)
  if argc > 1:
    if argv[1] == '-h':
      print(usage)
    for i in range(1, argc):
      examples.append(argv[i])

  else:
    for e in ['GaAs_fcc.py',
              'caffeine.py',
              'SnTe_2D.py',
              'CNO3.py',
              'example_qe_in.py',
              'example_poscar.py',
              'example_cif.py',
              'example_relax.py',
              'example_motion.py',
              'example_lammps_md.py',
              'example_cp2k.py',
              'example_charge_density.py',
              'example_BZ.py',
              'example_isosurface1.py',
              'example_isosurface2.py',
              'example_isosurface_clip.py',
              'example_colored_isosurface.py',
              'example_textured_isosurface.py',
              'example_textured_isosurface_clip.py']:
      examples.append(e)

  for e in examples:
    command = 'python {}'.format(e)
    print(command)
    if system(command) != 0:
      print('Aborting'); break
