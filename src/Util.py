import numpy as np

# Need to test this routine
def read_bxsf ( fname ):
  '''
  Read BXSF file, originally formatted for XCrysDen
  '''
  bands = None
  b_vec = None
  with open(fname) as f:
    l = f.readline()
    while l != '':
      if 'BANDGRID_3D_BANDS' in l:
        l = f.readline()
        nbnd = int(l.split()[0])
        nx,ny,nz = (int(v) for v in f.readline().split()); f.readline()
        b_vec = np.array([[float(v) for v in f.readline().split()] for _ in range(3)])
        bands = np.zeros((nx,ny,nz,nbnd), dtype=float)
        for n in range(nbnd):
          f.readline()
          for i in range(nx):
            for j in range(ny):
              ls = f.readline().split()
              for k in range(nz):
                bands[i,j,k,n] = float(ls[k])
        break
      l = f.readline()
  return b_vec,bands

# Need to test this routine
def read_xsf ( fname ):
  '''
  Read XSF file, originally formatted for XCrysDen
  '''
  pcoord = {}
  data = None
  nx = ny = nz = None
  pvec = np.empty((3,3), dtype=float)
  with open(fname) as f:
    l = f.readline()
    while l != '':
      if 'PRIMVEC' in l:
        for i in range(3):
          l = f.readline().split()
          for j in range(3):
            pvec[i,j] = float(l[j])
      elif 'PRIMCOORD' in l:
        l = f.readline().split()
        nat,nty = int(l[0]),int(l[1])
        for i in range(nat):
          l = f.readline().split()
          pos = [float(l[j]) for j in range(1,4)]
          if l[0] in pcoord:
            pcoord[l[0]].append(pos)
          else:
            pcoord[l[0]] = [pos]
      elif 'BEGIN_BLOCK_DATAGRID_3D' in l:
        l = f.readline()
        while not 'DATAGRID_3D' in l:
          l = f.readline()
        l = f.readline().split()
        nx,ny,nz = int(l[0]),int(l[1]),int(l[2])
        l = f.readline().split()
        origin = np.array([float(l[i]) for i in range(3)])
        for i in range(3):
          l = f.readline()
        ind = 0
        data = np.empty(nx*ny*nz, dtype=float, order='F')
        l = f.readline().split()
        while len(l) > 1:
          for v in l:
            data[ind] = float(v)
            ind += 1
          l = f.readline().split()
      l = f.readline()
  data = data.reshape((nx,ny,nz,1))
  return pvec,pcoord,data


def blender_xml ( fnpath, natoms, species, labels, atomic_poss, bonds, frame, cameras ):
  with open(fnpath, 'w') as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<XCP_DATA>\n')

    f.write('\t<SCENE>\n')
    if cameras is not None:
      for i,c in enumerate(cameras.values()):
        f.write('\t\t<CAMERA id="%d">\n'%i)
        f.write('\t\t\t<POSITION>%f %f %f</POSITION>'%tuple(c['position']))
        f.write('\t\t\t<ROTATION>%f %f %f</ROTATION>'%tuple(c['rotation']))
        f.write('\t\t</CAMERA>\n')
    f.write('\t</SCENE>\n')

    f.write('\t<SYSTEM>\n')
    for k,v in species.items():
      f.write('\t\t<SPECIES id="%d" label="%s">\n'%(species[k]['id'], k))
      f.write('\t\t\t<RADIUS>%f</RADIUS>\n'%species[k]['radius'])
      f.write('\t\t\t<COLOR>%f %f %f %f</COLOR>\n'%species[k]['color'])
      f.write('\t\t</SPECIES>\n')

    for i,a in enumerate(atomic_poss):
      key = labels[i%natoms]
      f.write('\t\t<ATOM id="%d" species="%d">'%(i,species[key]['id']))
      f.write('%f %f %f</ATOM>\n'%tuple(a))

    DEFAULT_BOND_TYPE = 1
    for i,b in enumerate(bonds):
      tup = (i, DEFAULT_BOND_TYPE, b[0], b[1])
      f.write('\t\t<BOND id="%d" type="%d" A="%d" B="%d"></BOND>\n'%tup)
    f.write('\t</SYSTEM>\n')

    vertices,edges = frame[0],frame[1]
    f.write('\t<FRAME>\n')
    for i,v in enumerate(vertices):
      f.write('\t\t<VERTEX id="%d">%f %f %f</VERTEX>\n'%((i,)+tuple(v)))
    for i,e in enumerate(edges):
      f.write('\t\t<EDGE id="%d" A="%d" B="%d"></EDGE>\n'%(i, e[0], e[1]))
    f.write('\t</FRAME>\n')

    f.write('</XCP_DATA>\n')
