from XtraCrysPy import XtraCrysPy as XCP

cpy = XCP.XtraCrysPy(inputfile='SnTe_data/SnTe.scf.in')
cpy.plot_spin_texture ( 'SnTe_data/Fermi_surf_band_4_0.npz', 'SnTe_data/spin_text_band_4.npz', colors=None, e_up=2, e_dw=-2, title='', w_width=1000, w_height=750, f_color=(1,1,1), bg_color=(0,0,0))
