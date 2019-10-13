from XCrysPy import XCrysPy as XCP

if '__main__' == __name__:

  # Spin texture plotted in reciprocal space
  # Arrow direction represents spin expectation and color represents energy.
  prefix = 'SnTe_data/'
  cpy = XCP.XCrysPy(qe_fname=prefix+'SnTe.scf.in')
  cpy.plot_spin_texture(prefix+'Fermi_surf_band_4_0.npz', prefix+'spin_text_band_5.npz')
