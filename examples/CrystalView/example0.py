from XtraCrysPy import XtraCrysPy,Model
import numpy as np

#struct = {'lattice':5.4*np.array([[-1,0,1],[0,1,1],[-1,1,0]]), 'species':['Si','Ge'], 'abc':np.array([[0.,0,0],[.25,.25,.25]]), 'colors':{'Si':[1,0,0], 'Ge':[1,1,1]}, 'bonds':{'Si_Ge':6}}
struct = {'lattice':5.4*np.array([[-1,0,1],[0,1,1],[-1,1,0]]), 'species':['H','U'], 'abc':np.array([[0.,0,0],[.25,.25,.25]]), 'colors':{'Si':[1,0,0]}, 'bonds':{'H_U':6}}
model = Model.Model(struct)

xcp = XtraCrysPy.XtraCrysPy()
xcp.render_atomic_model(model, nsc=(2,2,2), bond_type='primary')
xcp.start_crystal_view()

