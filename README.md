
# pysvr

A collection of python routines to calculate various parameters and indices relevant to severe convective storms. `pysvr` was initially created to fill a need to quickly calculate convective availible potential energy from 3-D/4-D gridded data using wrapped Fortran. 



## Acknowledgements
- Fortran CAPE code from RIP4 (NCL)
- SRH routines from [SHARPpy developers](https://github.com/sharppy/SHARPpy)
- Kim Hoogewind ([wx-pythonista](https://github.com/wx-pythonista))
## Authors

- [Victor Gensini](https://www.github.com/vgensini)


## Contributing

Contributions are always welcome!


## Example functions
(will hopefully build an API doc in the future...)

```python
import pysvr

#Calculate 100-hPa mixed layer CAPE
pysvr.thermodynamic.cape_plev_4D()

#Calculate the violent tornado parameter
pysvr.indicies.vtp()

#Calculate the bulk wind shear
pysvr.kinematic.bulk_wind_shear()
                                                                
```
