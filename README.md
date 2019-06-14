### csndm
```python
SSN.csndm(self, cells=None, normalize=1, to_csv=1, nodeW=0, *args, **kwargs)
```
Construction of network degree matrix
The function performs the transformation from gene expression matrix to network degree matrix (ndm).

Parameters:

    data    Gene expression matrix (TPM/RPKM/FPKM/count), rows = genes, columns = cells. otu_even.table
    alpha    Significant level (eg. 0.001, 0.01, 0.05 ...), Default = 0.01
    boxsize  Size of neighborhood, Default = 0.1 (nx(k) = ny(k) = 0.1*n)
    normalize \t1  result is normalized (Default);
               0  result is not normalized
Note:
    If gene expression matrix is sparse, use the sparse matrix will accelerate the calculation and reduce memory footprint
    data = sparse(data); upper = sparse(upper); lower = sparse(lower);
  
