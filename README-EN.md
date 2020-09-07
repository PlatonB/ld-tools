## Introduction.
_ld-tools_ is a self-written analogue of _LDlink_ and _HaploReg_, focused on non-server use. Created in _Python_ with _SQL_ elements.

| Program | Destination |
| ------- | ----------- |
| _ld_area_ | finds for each variant all variants inside two thresholds: not lower than a certain value of LD and not further than the specified flank size |
| _ld_triangle_ | builds matrices of LD values of all possible variant pairs |

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

## Restrictions.
- Programs skip multiallelic and non-rs variants.
- Working with X and Y chromosomes can be [problematic](https://github.com/samtools/bcftools/issues/1154).
- The first launch requires a stable Internet connection.
- Large _ld-triangle_'s heatmaps may not be rendered.
