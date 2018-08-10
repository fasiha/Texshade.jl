# Texshade

I downloaded the Version 5 CHGIS DEM (digital elevation elev) from https://doi.org/10.7910/DVN/M7WEFY, then
```bash
unzip v5_dem.zip
cd v5_dem
gdal_translate -of ENVI chgis_dem.tif chgis_dem.envi
cp chgis_dem.envi* ~/.julia/v0.6/Texshade
```
The metadata from `gdalinfo` is 
```
Size is 6299, 6170
Coordinate System is:
PROJCS["Xian_1980_GK_Zone_19",
    GEOGCS["GCS_Xian_1980",
        DATUM["Xian_1980",
            SPHEROID["Xian_1980",6378140.0,298.257]],
        PRIMEM["Greenwich",0.0],
        UNIT["Degree",0.0174532925199433]],
    PROJECTION["Transverse_Mercator"],
    PARAMETER["False_Easting",19500000.0],
    PARAMETER["False_Northing",0.0],
    PARAMETER["Central_Meridian",111.0],
    PARAMETER["Scale_Factor",1.0],
    PARAMETER["Latitude_Of_Origin",0.0],
    UNIT["Meter",1.0]]
Origin = (15781036.966107804328203,6567324.558187752030790)
Pixel Size = (1000.000000000000000,-1000.000000000000000)
Corner Coordinates:
Upper Left  (15781036.966, 6567324.558) ( 67d28'27.97"E, 48d 5'39.29"N)
Lower Left  (15781036.966,  397324.558) ( 79d18' 5.49"E,  3d 3'16.38"N)
Upper Right (22080036.966, 6567324.558) (149d40'36.28"E, 52d32'58.76"N)
Lower Right (22080036.966,  397324.558) (133d36'18.25"E,  3d18'57.26"N)
Center      (18930536.966, 3482324.558) (105d 1'19.30"E, 31d19'24.68"N)
Band 1 Block=64x64 Type=Int16, ColorInterp=Gray
  Min=-9983.000 Max=8596.000
  Minimum=-9983.000, Maximum=8596.000, Mean=-131.535, StdDev=2472.275
  Overviews: 3150x3085, 1575x1543, 788x772, 394x386, 197x193
Metadata:
  LAYER_TYPE=athematic
  STATISTICS_EXCLUDEDVALUES=-32768
  STATISTICS_MAXIMUM=8596
  STATISTICS_MEAN=-131.5354395341
  STATISTICS_MEDIAN=199
  STATISTICS_MINIMUM=-9983
```

Load the data into Julia
```julia

function readbin(filename, T)
    bytes = stat(filename).size
    fid = open(filename, "r")
    data = read(fid, T, bytes ÷ sizeof(T))
    close(fid)
    data
end
elev = convert(Array{Float32}, reshape(readbin("chgis_dem.envi", Int16), 6299, :))'

nodata = -32768f0
elev[elev .== nodata] = 0
```

```julia
import DSP
function texshadeBasic(elev, α=0.5; circular=false)
  pow = α * 0.5 # encode the sqrt here
  spectrum = rfft(elev)
  fx2s = DSP.fftfreq(size(elev, 2)).^2
  fy2s = DSP.rfftfreq(size(elev, 1)).^2
  maxAllowed = circular ? 0.5 : 1.0
  for col in 1:size(spectrum, 2)
    fx2 = fx2s[col]
    for row in 1:size(spectrum, 1)
      scale = min(maxAllowed, (fx2 + fy2s[row])^pow)
      spectrum[row, col] *= scale
    end
  end
  irfft(spectrum, size(elev, 1))
end

tex = texshadeBasic(elev)
texCirc = texshadeBasic(elev; circular=true)

import PyPlot
plt = PyPlot
plt.figure(); plt.imshow(elev); plt.colorbar()

plt.figure(); plt.imshow(tex); plt.colorbar(); plt.gci()[:set_clim](vec(quantile(vec(tex), [0.01 .99])))

plt.figure(); plt.imshow(texCirc); plt.colorbar(); plt.gci()[:set_clim](vec(quantile(vec(texCirc), [0.01 .99])))

plt.figure(); plt.imshow(tex - texCirc); plt.colorbar(); plt.title("tex - texCirc")

plt.figure(); plt.plt[:hist](vec(tex - texCirc), 128)


texQ = quantile(vec(tex), [0.01 .99])
texCircQ = quantile(vec(texCirc), [0.01 .99])

import PyPlot
function arrayToGrayPNG(arr, fname; lo=-Inf, hi=Inf)
  if (lo > -Inf)
    arr = max.(lo, arr)
  else
    lo = minimum(arr)
  end
  if (hi < Inf)
    arr = min.(hi, arr)
  else
    hi = maximum(arr)
  end
  
  output = convert.(UInt8, floor.(((arr - lo) / (hi - lo)) * 255.99))
  PyPlot.imsave(fname, output)
end

arrayToGrayPNG(tex, "tex.png"; lo=texQ[1], hi=texQ[2])
arrayToGrayPNG(texCirc, "texCirc.png"; lo=texCircQ[1], hi=texCircQ[2])


```

