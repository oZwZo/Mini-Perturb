# Mini-Perturb

## Clone repository and install prereqruisite
```
git clone https://github.com/oZwZo/Mini-Perturb.git
conda create -n GLM
conda activate GLM
conda install python>3.9
bash prerequisite.sh
```

## Download HoxB8 Data
```
bash data_download.sh
```

## Build GLM using Mini-Pert

```python
import os, sys
import scanpy as sc
from MiniPert.model import Perturb_NBGLM

# load adata
adata = sc.read_h5ad("<PATH-of-H5AD>")

# initiate model
# replace <> with your own key
# add interaction by passing variable pairs to `interaction_terms` 
pert_glm = Perturb_NBGLM(
    adata = adata, 
    perturb_key = "<PERTRUB>",
    plate_key = "<BATCH>",
    interaction_terms = None
)
pert_glm.fit(n_jobs=20)
```
