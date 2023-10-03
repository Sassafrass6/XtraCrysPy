from XtraCrysPy.Atomic import Atomic

a_info = {'radii':.1,
          'colors':{1:(1,0,0),2:(1,1,1)}}

xcp = Atomic(params=a_info, model='data_files/traj.20k.lmp', ftype='lammps-traj', md_perc=0.005)

# Animate with a duration of 50ms per frame, and no automatic restart
xcp.animate(fdt=50, restart_delay=-1)
xcp.start_crystal_view()
