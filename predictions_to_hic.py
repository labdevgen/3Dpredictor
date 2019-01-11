import pandas as pd
from matrix_plotter import MatrixPlotter
from matplot2hic import MatPlot2HiC

out_dir = "/mnt/scratch/ws/psbelokopytova/201901151331psbelokopytova/3DPredictor/out/hic_files"
fname = 
predicted_data = pd.read_csv("", sep="\t")
control_data=
mp = MatrixPlotter()
mp.set_data(predicted_data)
mp.set_control()
MatPlot2HiC(mp,fname, out_dir)
