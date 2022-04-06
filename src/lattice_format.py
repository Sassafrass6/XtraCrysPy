
def lattice_format_QE ( ibrav, celldm ):
  '''
    Given a space group and lattice, create a lattice in terms of QE's ibrav and compute atomic positions in this lattice

    Arguments:
      ibrav (int): ibrav index for pw.x
      celldm (list): Lattice parameters (as defined by quantum espresso)

    Returns:
      (ndarray): Lattice vectors
  '''
  import numpy as np

  a = celldm[0]
  boa = celldm[1]
  coa = celldm[2]
  cAB = celldm[3]
  cAC = celldm[4]
  cBC = celldm[5]

  nlat = None
  if ibrav == 0:
    raise ValueError('Can not format lattice for ibrav==0.')
  elif ibrav == 1:
    nlat = a * np.array([[1,0,0],[0,1,0],[0,0,1]])
  elif ibrav == 2:
    nlat = a/2 * np.array([[-1,0,1],[0,1,1],[-1,1,0]])
  elif ibrav == 3:
    nlat = a/2 * np.array([[1,1,1],[-1,1,1],[-1,-1,1]])
  elif ibrav == 4:
    nlat = a * np.array([[1,0,0],[-1/2,np.sqrt(3)/2,0],[0,0,coa]])
  elif abs(ibrav) == 5:
    tx = np.sqrt((1-cAB)/2)
    ty = np.sqrt((1-cAB)/6)
    tz = np.sqrt((1+2*cAB)/3)
    if ibrav == 5:
      nlat = a * np.array([[tx,-ty,tz],[0,2*ty,tz],[-tx,-ty,tz]])
    else:
      u = tz - 2*np.sqrt(2)*ty
      v = tz + np.sqrt(t)*ty
      nlat = a * np.array([[u,v,v],[v,u,v],[v,v,u]])
  elif ibrav == 6:
    nlat = a * np.array([[1,0,0],[0,1,0],[0,0,coa]])
  elif ibrav == 7:
    nlat = a/2 * np.array([[1,-1,coa],[1,1,coa],[-1,-1,coa]])
  elif ibrav == 8:
    nlat = a * np.array([[1,0,0],[0,boa,0],[0,0,coa]])
  elif abs(ibrav) == 9:
    nlat = a * np.array([[1/2,boa/2,0],[-1/2,boa/2,0],[0,0,coa]])
    if ibrav == -9:
      nlat[0,1] *= -1
      nlat[1,0] *= -1
  elif ibrav == 10:
    nlat = a/2 * np.array([[1,0,coa],[1,boa,0],[0,boa,coa]])
  elif ibrav == 11:
    nlat = a/2 * np.array([[1,boa,coa],[-1,boa,coa],[-1,-boa,coa]])
  elif abs(ibrav) == 12:
    if ibrav == 12:
      sAB = np.sin(np.arccos(cAB))
      nlat = a * np.array([[1,0,0],[boa*cAB,boa*sAB,0],[0,0,coa]])
    else:
      sAC = np.sin(np.arccos(cAC))
      nlat = a * np.array([[1,0,0],[0,boa,0],[coa*cAC,0,coa*sAC]])
  elif ibrav == 13:
    sAB = np.sin(np.arccos(cAB))
    nlat = a * np.array([[1/2,0,-coa/2],[boa*cAB,boa*sAB,0],[1/2,0,coa/2]])
  elif ibrav == 14:
    sAB = np.sin(np.arccos(cAB))
    t2 = coa * (cBC-cAC*cAB) / sAB
    t3 = coa * np.sqrt(1 + 2*cBC*cAC*cAB - cBC**2 - cAC**2 - cAB**2) / sAB
    nlat = a * np.array([[1,0,0],[boa*cAB,boa*sAB,0],[coa*cAC,t2,t3]])
  elif ibrav == 91:
    nlat = a * np.array([[1,0,0],[0,boa/2,-coa/2],[0,boa/2,coa/2]])

  return nlat
