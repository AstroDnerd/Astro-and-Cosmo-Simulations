import sys
import glob
import os
import imageio
from natsort import natsorted
from pdf2image import convert_from_path
from PIL import Image  
import PIL  

example_name = sys.argv[1]
plot_name = sys.argv[2]
path = '/home/astronerd/Desktop/AREPO/libraries/arepo-master/run/examples/'+example_name+'/plots/'+plot_name+'*.pdf'
files = glob.glob(path)
files = natsorted(files)
image_io = []
make_gif_command = 'gifmaker -i '
for i in range(len(files)):
	images = convert_from_path(files[i])
	images = images[0].save('/home/astronerd/Desktop/AREPO/test_run_plots/temp/'+str(i)+'.jpg')
	image_io.append(imageio.imread('/home/astronerd/Desktop/AREPO/test_run_plots/temp/'+str(i)+'.jpg'))
	print("done!"+str(i))

    
imageio.mimsave('/home/astronerd/Desktop/AREPO/test_run_plots/'+example_name+'_'+plot_name+'.gif', image_io, duration=0.5)
delete_temp = 'rm /home/astronerd/Desktop/AREPO/test_run_plots/temp/*'
os.system(delete_temp)