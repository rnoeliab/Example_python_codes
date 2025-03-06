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

lat_new = np.arange(-90, 90 + 1.5, 1.5)
lon_new = np.arange(-180, 180 + 1.5, 1.5)

# Cria o grid 2D de coordenadas
lon2d, lat2d = np.meshgrid(lon_new, lat_new)

Nt = co2.sizes["time"]
Nz = co2.sizes["level"]
Nlat = len(lat_new)
Nlon = len(lon_new)

co2_out = np.full((Nt, Nz, Nlat, Nlon), np.nan, dtype=np.float32)
co2_out = np.full((Nt, Nz, Nlat, Nlon), np.nan, dtype=np.float32)
co2_out = np.full((Nt, Nz, Nlat, Nlon), np.nan, dtype=np.float32)

lon_old = co2.Mesh_face_lon.values
lat_old = co2.Mesh_face_lat.values

for it in range(Nt):
    for iz in range(Nz):
        # Pega valores do CO2 para este tempo e nível (um array 1D de tamanho nMeshface)
        co2_values = co2.isel(time=it, level=iz).values
        # Interpola usando o método "nearest" (poderia ser "linear" ou "cubic" também)
        co2_interp = griddata(
            (lon_old, lat_old),
            co2_values,
            (lon2d, lat2d),
            method="nearest"
        )
        co2_out[it, iz, :, :] = co2_interp

co2_regridded = xr.DataArray(
    data=co2_out,
    coords={
        "time": co2.time,
        "level": co2.level,
        "lat": lat_new,
        "lon": lon_new
    },
    dims=("time", "level", "lat", "lon"),
    name="CO2"
)
co2_regridded.attrs = co2.attrs
co2_regridded.attrs["regridded_from"] = "unstructured to 1.5°x1.5° lat-lon"
co2_regridded = co2_regridded.to_dataset(name="co2")
#co2_regridded.to_netcdf("CO2_1p5deg.nc")


```

Here, there are some examples for you to learn python 
