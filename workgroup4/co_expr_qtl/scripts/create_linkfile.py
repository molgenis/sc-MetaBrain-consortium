import gzip
import pandas as pd

fh = open("./reffiles/sample_id.txt","rt")
elems = []
for line in fh:
    elems.append(line.strip())
fh.close()

print(f"{len(elems)} sample ids")

fho = open('./reffiles/linkfile.txt','w')
# header = fh.readline()
seen = set()

df = pd.read_csv("/groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2024-08-29-CombineDatasets/2024-09-02-BryoisEtAl/CombineFiles/wg0_file_directories.tsv",sep="\t")
pool_dataset = df.groupby('Pool')['Dataset'].apply(list).to_dict()

print(f"{len(pool_dataset)} pool-dataset links")

wctr = 0
for elem in elems:
    pool,ind = elem.split(".")
    if ind not in seen:
        seen.add(ind)
        dataset = pool_dataset[pool][0]
        fho.write(ind+"\t"+elem+"\t"+dataset+"\n")
        wctr += 1
fho.close()

print(f"{wctr} unique samples written")