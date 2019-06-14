Help on module csn:

NAME
    csn

DESCRIPTION
    python version of csn algorithm
    https://github.com/wys8c764/CSN

CLASSES
    builtins.object
        SSN
    
    class SSN(builtins.object)
     |  SSN(data, outdir='./', log=None)
     |  
     |  Construction of cell-specific networks
     |  模型构建过程用所有的样品数据，后续预测用整合有的大表做dm转化但仅输出少量样品cells.list的network和degree matrix
     |  在dblur features水平做矩阵融合
     |  The function performs the transformation from gene expression matrix to cell-specific network (csn).
     |      :param data: Gene expression matrix, rows = genes, columns = cells
     |  
     |  Methods defined here:
     |  
     |  __init__(self, data, outdir='./', log=None)
     |      default values when initialize. set log file
     |  
     |  csndm(self, cells=None, normalize=1, to_csv=1, nodeW=0, *args, **kwargs)
     |      Construction of network degree matrix
     |      ==========================================
     |      The function performs the transformation from gene expression matrix to network degree matrix (ndm).
     |          :param data: Gene expression matrix (TPM/RPKM/FPKM/count), rows = genes, columns = cells. otu_even.table
     |          :param alpha: Significant level (eg. 0.001, 0.01, 0.05 ...), Default = 0.01
     |          :param boxsize: Size of neighborhood, Default = 0.1 (nx(k) = ny(k) = 0.1*n)
     |          :param normalize: 1  result is normalized (Default); 0  result is not normalized
     |      If gene expression matrix is sparse, use the sparse matrix will accelerate the calculation and reduce memory footprint
     |      data = sparse(data); upper = sparse(upper); lower = sparse(lower);
     |      可用于机器学习，样品分类预测等
     |      只输出指定 cells 的degree matrix ，不指定就输出所有cell的全部genen dm
     |  
     |  csnet(self, cells=None, alpha=0.01, boxsize=0.1, edgeW=0, to_csv=0, *args, **kwargs)
     |      Construct the CSN for sepecified cells
     |      :param cells: Construct the CSNs for all cells, set cells = None (Default) otherwise input cells.list
     |      :param alpha: Significant level (eg. 0.001, 0.01, 0.05 ...)
     |             larger alpha leads to more edges, Default = 0.01
     |      :param boxsize: Size of neighborhood, Default = 0.1
     |      :param edgeW: 1  edge is weighted (statistic pxy(x))
     |             0  edge is not weighted (Default)
     |      :param nodeW: 1  node is weighted (gene or otu abundance)
     |             0  node is not wieghted (Default)
     |      :param csn: Cell-specific network, the kth CSN is in csn{k}
     |           rows = genes, columns = genes
     |      Too many cells or genes may lead to out of memory.
     |      
     |      :returns: cells in list format 
     |      :raises keyError: raises no exception
     |  
     |  get_cells(self, cells)
     |      :returns: cells in list format 
     |      :raises keyError: raises no exception
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)

FILE
    /home/yangfenglong/sphinx-apidoc/bin/csn.py
