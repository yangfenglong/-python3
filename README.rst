
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
