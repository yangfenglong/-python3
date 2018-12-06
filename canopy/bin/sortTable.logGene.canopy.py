#!/NJPROJ2/MICRO/PROJ/yangfenglong/software/miniconda3/bin/python3
# -*- coding: utf-8 -*-
# @Author: Yangfenglong
# @Date: 2018-08-10

import argparse
import pandas as pd
from datetime import datetime
from scipy import stats
import os
import multiprocessing
import numpy as np

def create_dir(dir):
    if not os.path.exists(dir):
        assert not os.system('mkdir %s' % dir)

class Canopy:
    def __init__(self, abuData, lenData): #default values when initialize
        self.dataset = abuData # gene abundance table: Unigenes.readsNum.even.xls
        self.lenData = lenData # gene length table: Unigenes.readsNum.screening.CDS.fa.len.xls
        self.t1 = 0.6
        self.t2 = 0.9
        self.min_step_cor = 0.99 #The value is hardcoded to not confuse users
        self.max_num_canopy_walks = 3  #The value is hardcoded to not confuse users
        self.minSize = 5  #filter the canopies which size < 5

    def setThreshold(self, t1, t2, method, minSize, outdir, outfile, p, nonzero, log):
        if t2 >= t1: # set parameters values
            self.t1 = t1
            self.t2 = t2
            self.method = method
            self.minSize = int(minSize)
            self.outdir = outdir
            self.outfile = outfile if outfile else 'cagCluster_{}_{}_{}.xls'.format(self.t1,self.t2,self.method)       
            self.p=p
            self.nonzero=nonzero
            self.log=log if log else 'cagCluster_{}_{}_{}.log'.format(self.t1,self.t2,self.method)
        else:
            print('t2 needs to be larger than t1!')        

#    def __iter__(self, data):
#        for line in open(data):
#            yield line.strip().split("\t")

    def correlationCoefficient(self, vec1, vec2):
        """caculate correlation coefficient rvalue"""
        nonzeroidx=np.union1d(vec1.nonzero()[0],vec2.nonzero()[0])
        vec1=vec1[nonzeroidx]
        vec2=vec2[nonzeroidx]
        if self.method =="spearman":
            return stats.spearmanr(vec1,vec2)[0]
        elif self.method == "pearson":
            return stats.pearsonr(vec1,vec2)[0]

    def recruitCanopy(self, neighbours, bait): #按照新的center重新招募canopy
        recruit_neighbours = []
        for neighbour in neighbours:
            cor = self.correlationCoefficient(bait, neighbour)
            if cor > self.t2:
                recruit_neighbours.append(neighbour)# 若相关系数大于t2，则将该点归入canopy
        #if(len(recruit_neighbours)==0): #没有招募到近邻
        #    return [neighbours, center]
        #else:                           #包括初始只有1个neighbour的 canopy，自己和自己cor=1
        df = pd.DataFrame(recruit_neighbours)
        center = df.median()
        return [recruit_neighbours, bait] if len(center.drop_duplicates())==1 else [recruit_neighbours, center]

    def canopyWalk(self, canopy):
        """walk 让每个canopy的质心中心化"""
        df = pd.DataFrame(canopy[1])
        center0=df.median()
        if len(center0.drop_duplicates())==1: #all values are seam
            [neighbours1, center1] = self.recruitCanopy(canopy[1], canopy[0]) #取出neighbour圈>t2的成员
        else:
            [neighbours1, center1] = self.recruitCanopy(canopy[1], center0)

        if len(neighbours1)<self.minSize : # 没有招募到近邻或近邻太少
            return 0
        else:
            [neighbours2, center2] = self.recruitCanopy(neighbours1, center1)   
        if len(neighbours2)<self.minSize :
            cnps = [cnp.name for cnp in neighbours1]
            return cnps
        cor = self.correlationCoefficient(center1, center2)       
        num_canopy_jumps_local = 0 #Max number of times the canopy will walk. This is a stop criterion for canopy walk.
        while (cor <= self.min_step_cor) and (num_canopy_jumps_local <= self.max_num_canopy_walks ): 
            """如果没超过循环阈值就一直walk下去，直到质心稳定"""
            center1 = center2
            neighbours1 = neighbours2
            [neighbours2,center2] = self.recruitCanopy(neighbours1,center1)
            if len(neighbours2)<self.minSize :
                cnps = [cnp.name for cnp in neighbours1]
                return cnps
            cor = self.correlationCoefficient(center1, center2)
            num_canopy_jumps_local += 1  
        if (len(neighbours1) > len(neighbours2)):   # Now we know that c1 and c2 are close enough and we should choose the one that has more neighbours
            cnps = [cnp.name for cnp in neighbours1]
            return cnps
        else:
            cnps = [cnp.name for cnp in neighbours2]
            return cnps

    def create_canopy(self):
        """main process to create canopies"""
        log = open(os.path.join(self.outdir, self.log), 'a') 
        log.write("{}: start to create canopies\n".format(datetime.now()))
        log.flush()

        """readin dataset and sort dataset by sum of rows :ascending=False"""
        df = pd.read_csv(self.dataset, sep='\t', index_col=0, header=0)
        sortindex=(df
            .apply(lambda x: x.sum(), axis=1)
            .sort_values(ascending=False)
            .index.get_values())
        sortdf=df.loc[sortindex]

        """1. initialize canopy"""
        canopies = []  # store cluster results >t1  
        candidates=[]  # store genes which (t1< cor <t2) with canopies that is uncertain whether it's a member in canopy
        geneNum=1
        line = sortdf.iloc[0]
        sortdf.drop(sortdf.index[0],inplace=True)
        while line[line!=0].size < self.nonzero or len(line.drop_duplicates())==1:
            line = sortdf.iloc[0] #不满足条件再读一行
            sortdf.drop(sortdf.index[0],inplace=True)
            geneNum += 1
        else:
            ori_center = line
            ori_center_list = [ori_center]   # 初始化ori_center点的canopy类容器
            canopies.append([ori_center,ori_center_list])
        canopyNum=1

        """2. create canopies"""
        for row in sortdf.iterrows():
            geneNum += 1
            line=row[1]
            if (line[line!=0].size > self.nonzero and len(line.drop_duplicates())!=1):
                NewPoint = line   # new point to be clusterd in the canopies
                counter = 0
                brk=0
                for canopy in canopies:
                    cor = self.correlationCoefficient(NewPoint, canopy[0])  # 计算选取的中心点P到每个点之间的相关系数
                    if cor >= self.t1: 
                        canopy[1].append(NewPoint)  # 若相关系数大于t1，则将该点归入canopy
                        counter +=1
                        if cor >= self.t2:
                            brk=1
                            break #下一个new point判断
                if not brk:
                    candidates.append(NewPoint)
                    if counter == 0:  # 没有并入之前的canopies: 新建一个Canopy
                        new_center = candidates.pop()
                        new_center_list = [new_center]
                        for i in range(len(candidates)-1,-1,-1):
                            cor = self.correlationCoefficient(candidates[i],new_center)
                            if cor >= self.t1:
                                new_center_list.append(candidates[i])
                                if cor >= self.t2:
                                    candidates.pop(i)
                        canopies.append([new_center,new_center_list])
                        canopyNum += 1
                        log.write("{}: create {} canopies from {} genes\n".format(datetime.now(),canopyNum,geneNum))
                        log.flush()

        log.write("{}: finish create canopies, start to centrolize the canopies\n".format(datetime.now()))
        log.flush()

        """3. centrolize and filter the canopies which size < self"""
        pool = multiprocessing.Pool(processes=self.p)
        centrolizedCanopies = []
        for canopy in (canopy for canopy in canopies if len(canopy[1]) > self.minSize):
            tmp = pool.apply_async(func=self.canopyWalk,args=(canopy,))
            if tmp.get():
                centrolizedCanopies.append(tmp.get())
        pool.close()
        pool.join()
        
        """4. uniq the centrolizedCanopies and del the subsets"""
        centrolizedCanopies.sort(reverse=True)
        uniqCanopies=[]
        i=0
        uniqCanopies.append(centrolizedCanopies[i])
        while i < len(centrolizedCanopies)-1:
            Subset=0
            for uniqCanopie in uniqCanopies:
                if set(uniqCanopie)>=set(centrolizedCanopies[i+1]):
                    i+=1
                    Subset=1
                    break
            if not Subset:
                uniqCanopies.append(centrolizedCanopies[i+1])
                i+=1

        log.write("{}: finish centrolize, filter the small canopies with size < {} and del the subsets\nstart to print out oricag.cluster file\n".format(datetime.now(),self.minSize))
        log.flush()

        """5. print out the oriCag file"""
        lenDict = {key:value for key,value in (line.strip().split('\t') for line in open(self.lenData))}
       # lenDict = {key:value for key,value in self.__iter__(self.lenData)}
        clusterid = 1
        create_dir(self.outdir)
        output = os.path.join(self.outdir, self.outfile)
        with open(output,'w') as f:
            f.write("cluster_id\tGenenum\tToatalLength\tLongest_gene_length\tLongest_gene_id\tCluster_genes\n")
            for cnp_genes in uniqCanopies:
                cnp_genes.sort()

                lens = [int(lenDict[gene]) for gene in cnp_genes]
                Genenum = len(cnp_genes)
                ToatalLength,Longest_gene_length = str(sum(lens)),str(max(lens))
                Longest_gene_id = [gene for gene in cnp_genes if int(lenDict[gene]) == max(lens)][0] 
                tmp = [str(clusterid),str(Genenum),ToatalLength,Longest_gene_length,Longest_gene_id,",".join(cnp_genes)]     
                f.write("\t".join(tmp) + "\n")
                clusterid += 1
        log.write("{}: finish print out the oriCag file\n".format(datetime.now()))
        log.flush()
        return clusterid-1

if __name__ == '__main__':
    ###################### argparse ###########################
    parser = argparse.ArgumentParser(
        description="canopy algorithm for co-abundance clustering of large gene catalogs cag analysis",
        usage="%(prog)s --help or -h",
        epilog="""
        Example:
           python  %(prog)s --abu Unigenes.readsNum.even.xls  --len Unigenes.readsNum.screening.CDS.fa.len.xls -m pearson --t1 0.6 --t2 0.9 --minSize 5 -p 10 --outdir ./ --outfile cag.ori.cluster.xls
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('--abu', '-a', required=True, help='Unigenes.readsNum.even.xls file (required)') 
    parser.add_argument('--len', '-l', required=True, help='Unigenes.readsNum.screening.CDS.fa.len.xls (required)')
    parser.add_argument('--method', '-m', default='pearson', choices=['spearman','pearson'], help='method for caculate correlation coefficient r_value, default = \'pearson\'') 
    parser.add_argument('--t1', default=0.6, type=float, help='min correlation r_value, default= 0.6')
    parser.add_argument('--t2', default=0.9, type=float, help='max correlation r_value, default= 0.9') 
    parser.add_argument('--minSize', default=5, help='filter the canopies which size <= minSize, default=5')
    parser.add_argument('-p', default=10, type=int, help='the process number for multiprocessing in canopy walk, default= 10')
    parser.add_argument('--nonzero', default=6, type=int, help='the number of samples whih contain the gene, for the min threthold when cadulate pearson correlation, default= 6')
    parser.add_argument('--outdir', '-d', default='./', help='a path to store analysis results, default=\'./\'')
    parser.add_argument('--outfile', '-o', help='the output file name')
    parser.add_argument('--logfile', help='the output log file name')
    args=parser.parse_args()
    ##########################################################

    t_s = datetime.now()    
    cnp = Canopy(args.abu, args.len) 
    cnp.setThreshold(args.t1, args.t2, args.method, args.minSize, args.outdir, args.outfile, args.nonzero, args.p, args.logfile)
    clusterNum = cnp.create_canopy()
    t_e = datetime.now()

    usedtime = t_e - t_s
    print('Get {} cag clusters, time used: [{}]'.format(clusterNum, usedtime))
