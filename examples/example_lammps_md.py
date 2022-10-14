from XtraCrysPy.Atomic import Atomic

a_info = {'radii':.1,
          'colors':{1:(1,0,0),2:(1,1,1)}}
xcp = Atomic(params=a_info, model='data_files/traj.20k.lmp', md_perc=0.005)
xcp.start_crystal_view()
