import numpy as np
from VertexTables import edge_offsets, corner_offsets, vertex_table

def get_cube_triangles ( index, grid, iso_val ):
  '''
    Evaluates the requred triangles within a single cube.

    Arugments:
      index (tuple of ints): A tuple with 3 indices for the current position in grid
      grid (ndarray): Array with the surface data to compare to iso_val
      iso_val (float): Value of the desired iso-surface

    Returns:
      (ndarray): List of 3d points, always in multiples of 3 points to define triangles.
  '''
  case = 0
  for i in range(8):
    ind = tuple(index[v]+corner_offsets[i][v] for v in range(3))
    grid_val = grid[ind]
    if grid_val - iso_val > 0:
      case |= 1<<i

  if case == 0 or case == 255:
    return

  vertex = 0
  vertices = []
  for i in range(5):
    for j in range(3):

      edge = vertex_table[case, vertex]
      if edge == -1:
        return np.array(vertices)

      v1 = np.array(index) + edge_offsets[edge, 0]
      v2 = np.array(index) + edge_offsets[edge, 1]

      s1 = grid[tuple(v1)] - iso_val
      s2 = grid[tuple(v2)] - iso_val

      diff = s1 - s2
      diff = .5 if diff==0. else s1/diff

      corner_pos = v1 + diff * (v2 - v1)

      vertices.append(corner_pos)

      vertex += 1


def marching_cubes ( grid, iso_val ):
  '''
    March through the entire grid and create a list of all triangles to draw

    Arguments:
      grid (ndarray): Array with the surface data to compare to iso_val
      iso_val (float): Value of the desired iso-surface

    Returns:
      (ndarray): List of 3d points representing ALL triangles that need to be drawn
  '''
  if len(grid.shape) != 3:
    print('Grid must be 3 dimensional.')
    quit()

  triangles = []
  n1,n2,n3 = grid.shape
  for i in range(n1-1):
    for j in range(n2-1):
      for k in range(n3-1):
        cube_tris = get_cube_triangles((i,j,k), grid, iso_val)
        if cube_tris is not None:
          if cube_tris.shape[0] > 3:
            for c in range(cube_tris.shape[0]//3):
              triangles.append(cube_tris[3*c:3*(c+1)])
          else:
            triangles.append(cube_tris)

  return np.array(triangles)
