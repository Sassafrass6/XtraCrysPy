from XtraCrysPy import XtraCrysPy,Model
import numpy as np

struct = {'lattice':np.array([[10,0,0],[0,10,0],[0,0,10.]]), 'species':['Si','Si'], 'abc':np.array([[0.,0,0],[.25,.25,.25]]), 'colors':{'Si':[1.,0,0]}, 'bonds':{'Si_Si':6}}
model = Model.Model(struct)

xcp = XtraCrysPy.XtraCrysPy()
xcp.render_atomic_model(model)
xcp.start_crystal_view()

