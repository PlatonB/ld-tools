## Introduction.
_ld-tools_ is a self-written analogue of LDlink and Haploreg, focused on non-server use. Created in _Python_ with _SQL_ elements.

| Program | Destination |
| ------- | ----------- |
| ld_area | finds for each of the requested variants all variants inside two thresholds: not lower than a certain value of LD and not further than the specified flank size |
| ld_triangle | builds matrices of LD values of all possible variant pairs of the input set |

## Dependencies.
The recommended method of dependency resolving is using Conda with Bioconda channel connected.
```
conda config --add channels defaults
```
```
conda config --add channels bioconda
```
```
conda config --add channels conda-forge
```
```
conda install python=3.7 pysam=0.15.4 plotly
```
