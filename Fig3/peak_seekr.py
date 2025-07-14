# seekr with vM25 lncRNA 


import numpy as np
import pandas as pd
from seekr.fasta_reader import Reader
from seekr.kmer_counts import BasicCounter
from seekr.pearson import pearson as seekrPearson

# generate norm vectors
bkg_norm_4 = BasicCounter('vM25.lncRNA.can.500.nodup.fa', k=4)
bkg_norm_4.get_counts()
mean_path = 'mean_4mers.npy'
std_path = 'std_4mers.npy'
np.save(mean_path, bkg_norm_4.mean)
np.save(std_path, bkg_norm_4.std)

# seekr count
allxp_count = BasicCounter(infasta='Xist_all_200cutpeaks.fa', outfile='Xist_allpeaks_200cut_4mers.csv', k=4,
                           mean='mean_4mers.npy', std='std_4mers.npy',
                           log2='Log2.post',label = True, binary=False)
allxp_count.make_count_file()


# pearson correlation
sim = seekrPearson(allxp_count.counts, allxp_count.counts)

# add row and column names to sim
pheaders = Reader('Xist_all_200cutpeaks.fa').get_headers()
# remove > from headers
pheaders = [header[1:] for header in pheaders]
df = pd.DataFrame(data=sim, index=pheaders, columns=pheaders)
# save the dataframe to csv
df.to_csv('Xist_allpeaks200_pearson.csv')


allxn_count = BasicCounter(infasta='Xist_null_200cut.fa', outfile='Xist_null_200cut_4mers.csv', k=4,
                           mean='mean_4mers.npy', std='std_4mers.npy',
                           log2='Log2.post',label = True, binary=False)
allxn_count.make_count_file()

# pearson correlation
sim = seekrPearson(allxp_count.counts, allxn_count.counts)

# add row and column names to sim
nheaders = Reader('Xist_null_200cut.fa').get_headers()
# remove > from headers
nheaders = [header[1:] for header in nheaders]
df = pd.DataFrame(data=sim, index=pheaders, columns=nheaders)
# save the dataframe to csv
df.to_csv('Xist_peakvsnull_pearson.csv')


#################Airn######################
allxp_count = BasicCounter(infasta='Airn_all_200cutpeaks.fa', outfile='Airn_allpeaks_200cut_4mers.csv', k=4,
                           mean='mean_4mers.npy', std='std_4mers.npy',
                           log2='Log2.post',label = True, binary=False)
allxp_count.make_count_file()


# pearson correlation
sim = seekrPearson(allxp_count.counts, allxp_count.counts)

# add row and column names to sim
pheaders = Reader('Airn_all_200cutpeaks.fa').get_headers()
# remove > from headers
pheaders = [header[1:] for header in pheaders]
df = pd.DataFrame(data=sim, index=pheaders, columns=pheaders)
# save the dataframe to csv
df.to_csv('Airn_allpeaks200_pearson.csv')

allxn_count = BasicCounter(infasta='Airn_null_200cut.fa', outfile='Airn_null_200cut_4mers.csv', k=4,
                           mean='mean_4mers.npy', std='std_4mers.npy',
                           log2='Log2.post',label = True, binary=False)
allxn_count.make_count_file()

# pearson correlation
sim = seekrPearson(allxp_count.counts, allxn_count.counts)

# add row and column names to sim
nheaders = Reader('Airn_null_200cut.fa').get_headers()
# remove > from headers
nheaders = [header[1:] for header in nheaders]
df = pd.DataFrame(data=sim, index=pheaders, columns=nheaders)
# save the dataframe to csv
df.to_csv('Airn_peakvsnull_pearson.csv')


#################Kot1######################
allxp_count = BasicCounter(infasta='Kot1_all_200cutpeaks.fa', outfile='Kot1_allpeaks_200cut_4mers.csv', k=4,
                           mean='mean_4mers.npy', std='std_4mers.npy',
                           log2='Log2.post',label = True, binary=False)
allxp_count.make_count_file()

# pearson correlation
sim = seekrPearson(allxp_count.counts, allxp_count.counts)

# add row and column names to sim
pheaders = Reader('Kot1_all_200cutpeaks.fa').get_headers()
# remove > from headers
pheaders = [header[1:] for header in pheaders]
df = pd.DataFrame(data=sim, index=pheaders, columns=pheaders)
# save the dataframe to csv
df.to_csv('Kot1_allpeaks200_pearson.csv')


allxn_count = BasicCounter(infasta='Kot1_null_200cut.fa', outfile='Kot1_null_200cut_4mers.csv', k=4,
                           mean='mean_4mers.npy', std='std_4mers.npy',
                           log2='Log2.post',label = True, binary=False)
allxn_count.make_count_file()

# pearson correlation
sim = seekrPearson(allxp_count.counts, allxn_count.counts)

# add row and column names to sim
nheaders = Reader('Kot1_null_200cut.fa').get_headers()
# remove > from headers
nheaders = [header[1:] for header in nheaders]
df = pd.DataFrame(data=sim, index=pheaders, columns=nheaders)
# save the dataframe to csv
df.to_csv('Kot1_peakvsnull_pearson.csv')
