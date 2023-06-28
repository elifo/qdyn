import pandas as pd
import os


CHUNKSIZE = int(1e8)
print ('CHUNKSIZE: ', CHUNKSIZE)
ICHUNK_STARTING_FROM = 10

quants_ox = ("t", "x", "y", "z", "v", "theta", "tau", "tau_dot", "slip", "sigma")
reader = pd.read_csv("output_ox", chunksize=CHUNKSIZE, header=None,
            names=quants_ox, delim_whitespace=True, comment="#")

for i_chunk, chunk in enumerate(reader):
    print ('i_chunk, chunk.shape', i_chunk, chunk.shape)
    if i_chunk > ICHUNK_STARTING_FROM:
       chunk.to_csv('./output_ox_chunk_'+str(i_chunk), sep=' ', header=None)
##
# only write the last chunk
#chunk.to_csv('./output_ox_chunk_'+str(i_chunk), sep=' ', header=None)

print ('Done!')
print ('*')
###