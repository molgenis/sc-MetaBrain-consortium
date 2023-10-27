#!/usr/bin/env python3

info_field_keys = None
with open("../1000G_b37/all_phase3_filtered.pvar", "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        info_field_keys = [field.split("=")[0] for field in line.strip("\n").split("\t")[-1].split(";")]
        break
f.close()
print("Looking for the following info fields:")
print(info_field_keys)
print("")

tmp_info_field_keys = info_field_keys.copy()
for field in tmp_info_field_keys:
    if "_" in field:
        a, b = field.split("_")
        info_field_keys.append(b + "_" + a)
del tmp_info_field_keys


print("Parsing pvar")
fh = open("20220422_3202_phased_SNV_INDEL_SV_b38_original.pvar", "r")
fho = open("20220422_3202_phased_SNV_INDEL_SV_b38.pvar", "w")
for i, line in enumerate(fh):
    if i % 250000 == 0:
        print("\t{} lines processed".format(i))
    if line.startswith("#"):
        fho.write(line)
    else:
        row_fields = line.strip("\n").split("\t")
        info_items = [item for item in row_fields[-1].split(";") if item.split("=")[0] in info_field_keys]
        fho.write("\t".join(row_fields[:-1]) + "\t" + ";".join(info_items) + "\n")
fh.close()
fho.close()