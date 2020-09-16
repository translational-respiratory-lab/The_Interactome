import pandas
import sys
import matplotlib.pyplot as plt
import seaborn as sns

df = pandas.read_csv(sys.argv[1],index_col=0)
df.columns=["Time to next exacerbation"]

sns.set(font_scale=0.1)
ax=sns.heatmap(df,linewidths=0.1, linecolor='black',cbar=False,center=0,vmin=0, vmax=1,cmap="seismic",square=True,xticklabels=False)
bottom, top = ax.get_ylim()
ax.set_ylim(bottom + 0.5, top - 0.5)
plt.savefig(sys.argv[2])
