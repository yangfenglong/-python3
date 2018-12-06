#help: /TJPROJ1/MICRO/yangfenglong/software/conda/miniconda3/bin/python3 heatmap.py kmeans9clusters_Alas2filt.csv  M1heatmap_Alas2filtdata.xls S1.heatmap.svg
from matplotlib import ticker, pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib as mpl
import pandas as pd #read tables
import sys
plt.switch_backend('agg') # does not plot to screen for RuntimeError: Invalid DISPLAY variable

# readin and sort tables
anno=pd.read_csv(sys.argv[1],sep=None,header=0,engine='python')
anno=anno.sort_values(by=[anno.columns.values[1],anno.columns.values[0]])
anno.set_index(anno.columns.values[0],inplace=True)
table=pd.read_table(sys.argv[2],sep=None,header=0,index_col=0,engine='python')
table=table[anno.index]

#print out the orderd tables
anno.to_csv("anno.xls",sep="\t")
table.to_csv("table.xls",sep="\t")

#for colors
rank={}
j=1
for i in anno.iloc[:,0].values:
    if i in rank:
        next
    else:
        rank[i]=j
        j +=1
anno.rank = [rank[i] for i in anno.iloc[:,0]]
#print(anno.rank)
# cmap=mpl.cm.Vega10
# colors set 
#colors = ['red', 'green', 'orange', 'blue', 'sandybrown', 'purple','pink','olive','cyan']
#https://matplotlib.org/examples/color/colormaps_reference.html  choose other cmps if you like
cmap = mpl.cm.tab20
colors = [mpl.colors.rgb2hex(cmap.colors[i-1])  for i in np.unique(anno.rank) if i <12]
cmap1 = mpl.colors.ListedColormap(colors)
#print(colors)
#for legend colors
col_legend = [colors[i-1]  for i in np.unique(anno.rank)]
label=[np.unique(anno.iloc[:,0])[i-1] for i in np.unique(anno.rank)]
#print(col_legend)
#print(label)
#for bar colors
color=[colors[i-1] for i in anno.rank]
#print(color)
# plot
fig,axes = plt.subplots(len(table)+1,1, sharex = True, sharey = False, figsize=(18,12),facecolor="white")
plt.subplots_adjust(left=None,bottom=None, right=0.8, top=None, wspace=0, hspace=0) 
# plot bars
for i in range(len(table)):
    ax=axes[i+1]
    table.iloc[i].plot(kind="bar",ax=ax,color=[color],width=1)
    ax.yaxis.set_label_position("right") 
    ax.set_ylabel(ylabel=table.iloc[i].name,rotation="horizontal",horizontalalignment="left")
    ax.set_autoscaley_on(True)
    ax.set_xticks([])
    ax.spines["bottom"].set_color('gray')
    ax.spines["top"].set_color('gray')
    ax.spines["left"].set_color('gray')
    ax.spines["right"].set_color('gray')   
    ax.set_yticks( np.arange(4) )
    ax.axes.title.set_visible(False)
    ax.tick_params(bottom=False,top=False,left=True,right=True,
                   labelbottom=False,labelleft=False,labelright=False,which="both")
    ax.tick_params(direction="in")

# top bar plot
ax=axes[0]
#ax.axes.get_xbound()
#ax.axes.get_ybound()[1]/2
extent = (ax.axes.get_xbound()[0],ax.axes.get_xbound()[1],ax.axes.get_ybound()[0],ax.axes.get_xbound()[1]*0.01)
ax.imshow([anno.rank], interpolation='none', cmap=cmap1, alpha=1, aspect="equal", extent=extent)#,norm=norm
ax.set_axis_off()

# for legend
patches =[mpatches.Patch(color=col_legend[i-1],label=np.unique(anno.iloc[:,0])[i-1]) for i in np.unique(anno.rank)]
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width , box.height])
ax.legend(handles=patches, bbox_to_anchor=(1.16, 0.2),ncol=1,borderaxespad=0,fontsize=10)
fig.savefig(fname=sys.argv[3])
