# Simple python examples

## Create a environment.yml
```
> conda activate py39
> conda env export --no-builds> environment.yml
```

Use this environment to run these examples:

```
a. Have Anaconda installed ([Installing Anaconda](https://github.com/rnoeliab/Installing_anaconda))
b. conda env create -f environment.yml
c. conda activate py39 
```

## Read CO2 data from CAMS-inversion
Since year 2019, CAMS has also been producing atmospheric inversions based on satellite retrievals from the
second Orbiting Carbon Observatory (OCO-2, https://oco.jpl.nasa.gov/oco-2-data-center/). Only retrieval data
over land are used.

### How to read the CO2 data from CAMS-inversion in python
```
import xarray as xr
import numpy as np
from scipy.interpolate import griddata

files = 'path_file'
ds = xr.open_dataset(files + 'cams73_latest_co2_conc_satellite_inst_202301.nc')
print(ds)
```
![co2](https://github.com/rnoeliab/Example_python_codes/blob/main/figures/co2_dataset.png)

```
ds.Mesh_face_lon.data.min(), ds.Mesh_face_lon.data.max()
co2 = ds.CO2.copy()
height_are = ds.height_above_reference_ellipsoid.copy()
Psurf = ds.Psurf.copy()
ap = ds.ap.copy()
bp = ds.bp.copy()
```

Here, there are some examples for you to learn python 
