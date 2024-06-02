
### Overview
GeoJikuu is an open-source Python library designed for geospatial analysis. Although it works in any Python environment, it is particularly suited for use with Pandas DataFrames in Jupyter Notebooks. The current version of GeoJikuu is considered alpha, with GeoJikuu v1.0 set for release later in 2024. Despite being an alpha version, it already offers functionalities for calculating spatial and spatiotemporal descriptive statistics, performing spatial and spatiotemporal aggregation, and conducting spatial and spatiotemporal hypothesis testing. Preliminary documentation is available at [geojikuu.com/docs](geojikuu.com/docs).

### Installation
GeoJikuu can be installed via: 
```python
pip install geojikuu
```

### Usage Example
Projecting coordinates to Cartesian form and running Global Moran's I:
```python
from geojikuu.hypothesis_testing.autocorrelation import GlobalMoranI
from geojikuu.preprocessing.projection import CartesianProjector

cartesian_projector = CartesianProjector("wgs84")

data = {
    "lat": [34.6870676, 34.696109, 34.6525807, 35.7146509, 35.6653623, 35.6856905, 
            33.5597115, 33.5716997, 33.5244701, 33.5153417, 33.5206116, 33.4866878],
    "lon": [135.5237618, 135.5121774, 135.5059984, 139.7963897, 139.7254906, 139.7514867,
            130.3818748, 130.4030704, 130.4063441, 130.4373212, 130.4841434, 130.5220605],
    "value": [2, 3, 1, 4, 4, 2, 5, 6, 5, 7, 8, 8]
}

df = pd.DataFrame.from_dict(data)

results = cartesian_projector.project(list(zip(df["lat"], df["lon"])))
df["cartesian_coordinates"] = results["cartesian_coordinates"]
unit_conversion = results["unit_conversion"]

global_moran_i = GlobalMoranI(data=df, coordinate_label="cartesian")
global_moran_i.run(input_field="value", critical_distance=10/unit_conversion)
```

### Future Updates
The following additions and improvements can be expected in GeoJikuu v1.0:
 - In-built visualisation
 - Optimised algorithms for efficient operations on large datasets
 - Many additional modules for performing various types of spatial, temporal, and spatiotemporal analyses
 - Improved packaging and module designs that are more suitable for representing and working with spatial structures
 - Improved documentation and the introduction of a blog and newsletter

### Contact
If you are interested in the project and have any questions or implementation requests, please direct your correspondence to admin@gaiaabstract.com.
Feedback from academic and professional researchers is particularly welcome. 

