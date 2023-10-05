import argparse

def parse_arguments ( ):

  parser = argparse.ArgumentParser('XtraCrysPy: A visualization tool for the inputs and outputs of electronic structure and molecular dynamics softwares.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-v', '--verbose', action='count', default=0, help='Set the verbosity (0 for low, 1 for high)')

  fn_group = parser.add_mutually_exclusive_group()
  fn_group.add_argument('-fn', '--filename', help='Path to file containing the atomic system')

  parser.add_argument('-ff', '--fileformat', help='Format of the provided filename (ASE standard)')

  parser.add_argument('-mf', '--multiframe', default=False, help='Multiple frames should be read from the provided filename')

  parser.add_argument('-bg', '--background', default=(0,0,0), help='Tuple containing (r,g,b) color for the background')

  parser.add_argument('-bl', '--bondlength', help='Maximum length for drawing bonds between atomic sites')
  return parser.parse_args()


def process_arguments ( args ):

  # Parse MultiFrame boolean argument
  if type(args.multiframe) == str and args.multiframe[0].isdigit():
    args.multiframe = int(args.multiframe)
    args.multiframe = bool(args.multiframe)

  # Parse Background Color
  if type(args.background) == str:
    iscap = lambda s : False if s in ['(',')','[',']'] else True
    args.background = tuple(float(v) for v in ''.join(filter(iscap,args.background)).split(','))
    for c in args.background:
      if c > 1 or c < 0:
        raise ValueError('Background r,b,g colors should be in the range [0,1]')

  # Parse bond length
  if args.bondlength is not None:
    args.bondlength = float(args.bondlength)

  return args

def xcp_main ( ):

  args = process_arguments(parse_arguments())

  params = {}
  if args.bondlength is not None:
    params['bonds'] = args.bondlength

  from .Atomic import Atomic
  xcp = Atomic(model=args.filename, background=args.background, params=params)
  xcp.start_crystal_view()
