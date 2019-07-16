"""
    Demo script to generate fractals and save them as NRRD file.

"""

from tlib.lung import *
import tlib.utils as tut
import nrrd

sf = SukiFractal(n=2.8,
                 r=.4,
                 d0=50,
                 max_generation=3,
                 volume_size=(1000, 1000, 1000),
    )

ix = 1
while True:
    file = 'lung_' + str(ix) + '.nrrd'
    if os.path.isfile(file) is True:
        ix += 1
    else:
        nrrd.write(file, sf.volume.voxels)
        break
#         tut.write_zipped_pickle(sf.volume.voxels, file)

print('Saved {}.'.format(file))
