import pandas
sim=pandas.read_csv("diff_75%-reduction.csv default edge.csv",index_col=0)
post=pandas.read_csv("diff_post.csv",index_col=0)

sim_nam=sim[["name","interaction"]].values
post_nam=post[["name","interaction"]].values

sim_nam=set([str(i) for i in sim_nam])
post_nam=set([str(i) for i in post_nam])

val=(len(post_nam.intersection(sim_nam)))/(len(sim_nam))
print(val*100)

val=(len(post_nam.intersection(sim_nam)))/(len(post_nam))
print(val*100)
