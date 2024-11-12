import gzip
import gzip
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, type=str,  help="")
parser.add_argument("--output", required=True, type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

fin = gzip.open(args.input,"rt")
fout = open(args.output,"wt")
samples = fin.readline().strip().split("\t")[1:]
for id in samples:
    fout.write(f"{id}\n")

fin.close()
fout.close()
print("Done.")