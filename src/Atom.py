from vpython import sphere,vector
class Atom:
  def __init__ ( self, pos, col=vector(1,1,1), species='' ):
    '''
    Initialize an Atom object

    Arguments:
      pos (vpython.vector): Atom position
      col (vpython.vector): Atom color
      speciec (str): Species identification string
    '''
    self.pos = pos
    self.col = col
    self.species = str(species)
    self.vpy_sph = sphere(pos=self.pos, color=self.col, radius=.75)
