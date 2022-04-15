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
    for e in ['silicon_fcc.py',
              'caffeine.py',
              'SnTe_2D.py',
              'example_qe_in.py',
              'example_poscar.py',
              'example_cif.py',
              'example_relax.py',
              'example_charge_density.py',
              'example_BZ.py',
              'example_fermi_surface1.py',
              'example_fermi_surface2.py',
              'example_fermi_surface_clip.py',
              'example_colored_fermi_surface.py',
              'example_textured_fermi_surface_clip.py']:
      examples.append(e)

  for e in examples:
    command = 'python {}'.format(e)
    print(command)
    if system(command) != 0:
      print('Aborting'); break
