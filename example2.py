from CrysPy import CrysPy

crystal = CrysPy(qe_fname='Fe.scf.in')
crystal.draw_cell(boundary=True)
