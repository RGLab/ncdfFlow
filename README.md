# ncdfFlow: A package that provides HDF5 based storage for cytometry data.

This package extends the `flowCore` infrastructure by storing the large volume of event-level data on disk as `HDF` format and only keeps the file handler and meta data in memory. Thus the memory consumption is significantly reduced.

### INSTALLATION




```r
# First, install it from bionconductor so that it will pull all the dependent packages automatically
library(BiocInstalller)
bicLite(ncdfFlow) 

# or install the latest version from github using devtools package 
install.packages("devtools") 
library(devtools) #load it
install_github("RGLab/ncdfFlow", ref="trunk")
```

### Unix/Linux/Mac users

To build the ncdfFlow package from source, make sure that HDF5 Library is present on your system:

If HDF5 is installed to some non-standard location, you may pass the environment variable --with-hdf5 to point to the correct location of HDF5, for example,
```bash
#install from github
install_github('RGLab/ncdfFlow', ref='trunk', args='--configure-args="--with-hdf5=<path-to-hdf>"')
 
#or install from locally downloaded tar ball
R CMD INSTALL ncdfFlow_x.y.z.tar.gz --configure-args="--with-hdf5='<path-to-hdf>'"
```
under '/path/to', you should find "include" and "lib" sub-folders that contain HDF5 headers and shared libraries. 

Also, make sure add the path of `libhdf5.so` (should be `lib` subfolder of `<path-to-hdf>`) to your environment variable `LD_LIBRARY_PATH` so that it can be found at runtime.

```bash
export LD_LIBRARY_PATH=<path-to-hdf>/lib:LD_LIBRARY_PATH
```

### Create `ncdfFlowSet` object


```r
library(ncdfFlow)

#read from FCS files
path <- system.file("extdata","compdata","data",package="flowCore")
files <- list.files(path,full.names=TRUE)[1:3]
fs <- read.ncdfFlowSet(files=files) #equivalent to flowCore::read.flowSet


#or convert the existing flowSet into ncdfFlowSet
data(GvHD)
fs <- GvHD[1:4]
fs <- ncdfFlowSet(fs)
fs
```

```
## An ncdfFlowSet with 4 samples.
## NCDF file : /tmp/RtmpENFQmY/ncfs85d52869372e.nc 
## An object of class 'AnnotatedDataFrame'
##   rowNames: s5a01 s5a02 s5a03 s5a04
##   varLabels: Patient Visit ... name (5 total)
##   varMetadata: labelDescription
## 
##   column names:
##     FSC-H, SSC-H, FL1-H, FL2-H, FL3-H, FL2-A, FL4-H, Time
```

### Use it as the same way as `flowSet` (except it is memory efficient and fast)

```r
pData(fs)
```

```
##       Patient Visit Days Grade  name
## s5a01       5     1   -6     3 s5a01
## s5a02       5     2    0     3 s5a02
## s5a03       5     3    6     3 s5a03
## s5a04       5     4   12     3 s5a04
```

```r
sampleNames(fs)
```

```
## [1] "s5a01" "s5a02" "s5a03" "s5a04"
```

```r
keyword(fs,"FILENAME")
```

```
##       FILENAME
## s5a01 "s5a01" 
## s5a02 "s5a02" 
## s5a03 "s5a03" 
## s5a04 "s5a04"
```

```r
colnames(fs)
```

```
## [1] "FSC-H" "SSC-H" "FL1-H" "FL2-H" "FL3-H" "FL2-A" "FL4-H" "Time"
```

```r
length(fs)
```

```
## [1] 4
```

```r
fs[[1]]
```

```
## flowFrame object 's5a01'
## with 3420 cells and 8 observables:
##      name              desc range minRange maxRange
## $P1 FSC-H        FSC-Height  1024        0     1023
## $P2 SSC-H        SSC-Height  1024        0     1023
## $P3 FL1-H         CD15 FITC  1024        1    10000
## $P4 FL2-H           CD45 PE  1024        1    10000
## $P5 FL3-H        CD14 PerCP  1024        1    10000
## $P6 FL2-A              <NA>  1024        0     1023
## $P7 FL4-H          CD33 APC  1024        1    10000
## $P8  Time Time (51.20 sec.)  1024        0     1023
## 153 keywords are stored in the 'description' slot
```

```r
fs[2:3]
```

```
## An ncdfFlowSet with 2 samples.
## NCDF file : /tmp/RtmpENFQmY/ncfs85d52869372e.nc 
## An object of class 'AnnotatedDataFrame'
##   rowNames: s5a02 s5a03
##   varLabels: Patient Visit ... name (5 total)
##   varMetadata: labelDescription
## 
##   column names:
##     FSC-H, SSC-H, FL1-H, FL2-H, FL3-H, FL2-A, FL4-H, Time
```


