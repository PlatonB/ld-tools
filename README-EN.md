## Introduction.
_ld-tools_ is a self-written analogue of _LDlink_ and _HaploReg_, focused on non-server use. Created in _Python_ with _SQL_ elements.

| Program | Main features |
| ------- | ------------- |
| _ld_area_ | selects for each variant all variants inside two thresholds: not lower than a certain value of LD and not further than the specified flank size |
| _ld_lite_ | outputs LD and distance for two variants |
| _ld_triangle_ | builds matrices of LD values of all possible variant pairs |

## Dependencies.
The recommended method of dependency resolving is using _Conda_ with _Bioconda_ channel connected.
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
conda install python=3.7 pysam=0.15.4 tabulate plotly
```

## Deployment.
1. [Download](https://github.com/PlatonB/ld-tools/archive/master.zip) the archive with the program.
2. Unpack it to any place.
3. In the terminal, navigate to the unzipped folder.
```
cd /path/to/ld-tools-master
```
4. Print the help for the component you are interested in.
```
python3 ld_triangle -h
```

## Restrictions.
- Programs skip multiallelic and non-rs variants.
- Working with X and Y chromosomes can be [problematic](https://github.com/samtools/bcftools/issues/1154).
- The first launch requires a stable Internet connection.
- Large _ld-triangle_'s heatmaps may not be rendered.

## Input.
- Folder with plain tables. They must contain a column with rsIDs. The position of this column can be any. The number of headers is arbitrary, but must be uniform between the tables.
- Folder for 1000 Genomes data. The data shall be downloaded and processed automatically at the first launch.

## Output.
### ld_area.
A separate file is created for each variant. Output file contains the source variant and all variants selected for it according to the specified conditions.
- TSV. Default choice. Contains r2, D' and the distance between found and requested variant, plus basic annotations of each variant.
- JSON. It is similar to TSV in content, but is more scripting oriented.
- Column of rsIDs. Minimalistic output for annotation from scratch.

### ld_triangle.
- LD matrix as heatmap. For advanced usage it is possible to derive diagram object in JSON format.
- LD matrix as table.
