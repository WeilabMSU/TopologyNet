# TopologyNet
Software for generating topological features for predictions of protein-ligand binding affinity and mutation induced protein stability changes.

## Protein-ligand binding
<img src="https://github.com/WeilabMSU/TopologyNet/blob/master/fig/binding_figure.PNG" width="500">

### Usage
See the folder ``Protein-ligand-binding``. The source codes are in the folder ``Protein-ligand-binding/TopBio``. An example protein-ligand complex structure for ``PDB:1A8I`` is provided. See the file ``test_example.py`` for generating the various topological descriptors for the example protein-ligand complex.

### Dependency
1. [Javaplex 4.3.1](https://github.com/appliedtopology/javaplex/releases/tag/4.3.1)
2. [R-TDA](https://cran.r-project.org/web/packages/TDA/index.html)
3. Numpy, Pickle

### Reference
[Cang, Zixuan, Lin Mu, and Guo-Wei Wei. "Representability of algebraic topology for biomolecules in machine learning based scoring and virtual screening." PLoS computational biology 14.1 (2018): e1005929.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005929)

[Cang, Zixuan, and Guo-Wei Wei. "TopologyNet: Topology based deep convolutional and multi-task neural networks for biomolecular property predictions." PLoS computational biology 13.7 (2017): e1005690.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005690)
