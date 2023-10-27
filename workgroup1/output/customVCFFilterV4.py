import gzip
import sys

if len(sys.argv) < 3:
    print("Usage: VCFin.vcf.gz vcfOutPrefix [samplelist]")
    sys.exit(0)

vcfIn = sys.argv[1]
vcfOut = sys.argv[2]

sampleFilterSet = None  # <-- for example to select EUR samples only
if len(sys.argv) > 3:
    sampleFilterSet = set()
    print("Reading samples to filter from: " + sys.argv[3])
    fh = open(sys.argv[3], 'r')
    for line in fh:
        sampleFilterSet.add(line.strip())
    fh.close()
    print("{} samples selected".format(len(sampleFilterSet)))

# variables
tresh_GQ = 20  # genotype quality
tresh_AB_lower = 0.2  # allelic balance lower cutoff
tresh_AB_upper = 0.8  # allelic balance upper cutoff
tresh_IB = -0.3  # inbreeding coefficient threshold
check_VQSR_SNV = True
tresh_VQSR_SNV = 99.8  # VQSR lower threshold for SNPs
check_VQSR_Indel = True
tresh_VQSR_Indel = 99.95  # VQSR lower bound for indels
remove_multiAllelic = True
remove_nonPASS_indel = True
remove_nonPASS_SNV = True
filter_lowcomplexity = True
tresh_MAF = 0.01  # minor allele frequency cutoff
tresh_CR = 0.99  # call rate threshold
tresh_HWE = 1E-6  # Hardy-Weinberg p-value cutoff
tresh_AD = 10  # minimum allelic depth
tresh_DP = 10  # minimum total depth

# output vars
stripINFOcol = True  # clear INFO column, because these values are likely invalid after filtering

debug = False
stopafterlines = 100000
printhwe = False
chrpos = None


def parseVQSR(vqsrstring, isIndel):
    if not vqsrstring.startswith("VQSR"):
        return True
    vqsrstring = vqsrstring.replace("VQSRTrancheINDEL", "")
    vqsrstring = vqsrstring.replace("VQSRTrancheSNP", "")
    vqsrstring = vqsrstring.replace("+", "")
    velems = vqsrstring.split("to")
    lower = float(velems[0])
    if isIndel:
        return (lower >= tresh_VQSR_Indel)
    else:
        return (lower >= tresh_VQSR_SNV)


def summarizeSettings():
    outln = "##" + sys.argv[0]
    outln += (" tresh_GQ =" + str(tresh_GQ) +
              ";tresh_AB_lower =" + str(tresh_AB_lower) +
              ";tresh_AB_upper =" + str(tresh_AB_upper) +
              ";tresh_IB =" + str(tresh_IB) +
              ";check_VQSR_SNV =" + str(check_VQSR_SNV) +
              ";tresh_VQSR_SNV =" + str(tresh_VQSR_SNV) +
              ";check_VQSR_Indel =" + str(check_VQSR_Indel) +
              ";tresh_VQSR_Indel =" + str(tresh_VQSR_Indel) +
              ";remove_multiAllelic =" + str(remove_multiAllelic) +
              ";remove_nonPASS_indel =" + str(remove_nonPASS_indel) +
              ";remove_nonPASS_SNV =" + str(remove_nonPASS_SNV) +
              ";filter_lowcomplexity =" + str(filter_lowcomplexity) +
              ";tresh_MAF =" + str(tresh_MAF) +
              ";tresh_HWE =" + str(tresh_HWE) +
              ";tresh_CR =" + str(tresh_CR) +
              ";tresh_AD =" + str(tresh_AD) +
              ";tresh_DP =" + str(tresh_DP))
    outln += "\n"
    return outln


def calculateHWE(obs_hets, obs_hom1, obs_hom2):
    # still need to implement this...
    obs_homc = obs_hom1
    if obs_hom1 < obs_hom2:
        obs_homc = obs_hom2
    obs_homr = obs_hom2
    if obs_hom1 < obs_hom2:
        obs_homr = obs_hom1
    rare_copies = 2 * obs_homr + obs_hets
    l_genotypes = obs_hets + obs_homc + obs_homr
    if printhwe:
        print("{}\t{}".format(rare_copies, l_genotypes))
    if l_genotypes == 0:
        return -1

    het_probs = [0] * (rare_copies + 1)
    mid = int(rare_copies) * int(2 * l_genotypes - rare_copies) / int(
        2 * l_genotypes)
    mid = int(mid)
    if printhwe:
        print("{}".format(mid))

    if mid % 2 != rare_copies % 2:
        mid += 1
    mid = int(mid)
    if printhwe:
        print("{}".format(mid))
    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = l_genotypes - curr_hets - curr_homr
    if printhwe:
        print("{}\t{}\t{}".format(curr_hets, curr_homr, curr_homc))
    het_probs[int(mid)] = 1.0
    sum = het_probs[int(mid)]

    curr_hets = int(mid)
    while curr_hets > 1:
        het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (
                    curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (
                    curr_homc + 1.0))
        sum += het_probs[curr_hets - 2]
        if printhwe:
            print("{}\t{}".format(curr_hets, het_probs[curr_hets - 2]))
        curr_homr += 1
        curr_homc += 1
        curr_hets -= 2

    curr_hets = int(mid)
    curr_homr = (rare_copies - mid) / 2
    curr_homc = l_genotypes - curr_hets - curr_homr
    if printhwe:
        print("{}\t{}\t{}".format(curr_hets, curr_homr, curr_homc))
    while curr_hets <= (rare_copies - 2):
        het_probs[curr_hets + 2] = het_probs[
                                       curr_hets] * 4.0 * curr_homr * curr_homc / (
                                               (curr_hets + 2.0) * (
                                                   curr_hets + 1.0))
        sum += het_probs[curr_hets + 2]
        if printhwe:
            print("{}\t{}".format(curr_hets, het_probs[curr_hets - 2]))
        curr_homr -= 1
        curr_homc -= 1
        curr_hets += 2

    i = 0
    while i <= rare_copies:
        het_probs[i] /= sum
        if printhwe:
            print("{}\t{}\t{}".format(i, het_probs[i], sum))
        i += 1

    p_hwe = 0.0
    i = 0
    while i <= rare_copies:
        if het_probs[i] <= het_probs[obs_hets]:
            p_hwe += het_probs[i]
        i += 1
    if p_hwe > 1:
        p_hwe = 1
    return p_hwe


def getMAF(dosages):
    called = 0
    nrhomA = 0
    nrhets = 0
    nrhomB = 0
    nrA = 0
    nrB = 0

    for dosage in dosages:
        if dosage > -1:
            called += 1
            if dosage == 0:
                nrhomA += 1
                nrA += 2
            elif dosage == 1:
                nrhets += 1
                nrA += 1
                nrB += 1
            else:
                nrhomB += 1
                nrB += 2
    if called == 0:
        return 0
    maf = nrA / (called * 2)
    if maf > 0.5:
        maf = 1 - maf
    return maf


def getStats(dosages, tresh_MAF, tresh_CR, tresh_HWE):
    called = 0
    nrhomA = 0
    nrhets = 0
    nrhomB = 0
    nrA = 0
    nrB = 0
    for dosage in dosages:
        if dosage > -1:
            called += 1
            if dosage == 0:
                nrhomA += 1
                nrA += 2
            elif dosage == 1:
                nrhets += 1
                nrA += 1
                nrB += 1
            else:
                nrhomB += 1
                nrB += 2
    varOk = True
    maf = 0
    hwe = -1
    callrate = 0
    if called == 0:
        varOk = False
    else:
        callrate = called / len(dosages)
        maf = nrA / (called * 2)
        if maf > 0.5:
            maf = 1 - maf
        if callrate >= tresh_CR and maf >= tresh_MAF:
            hwe = calculateHWE(nrhets, nrhomA,
                               nrhomB)  # only calculate HWE when MAF and CR above thresholds
            if printhwe:
                print("HWE for variant: " + chrpos + "\t" + str(hwe))
                print(
                    "nrhets: {}, nrhomA: {}, nrhomB: {}".format(nrhets, nrhomA,
                                                                nrhomB))
                sys.exit(0)
            if hwe < tresh_HWE:
                varOk = False
        else:
            varOk = False
    statout = "CR={};MAF={};HWE={:0.3e}".format(round(callrate, 3),
                                                round(maf, 3), hwe)
    return [varOk, statout]


def parseLine(lineNumber, line, includeSample):
    # line with genotypedata
    elems = line.strip().split("\t", 10)

    # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
    # 0     1       2       3       4       5       6       7       8
    chr = elems[0]
    pos = elems[1]
    global printhwe
    global chrpos
    chrpos = chr + ":" + pos
    #       if chrpos == "1:14590":
    #               printhwe = True
    id = elems[2]
    logid = chr + ":" + pos + ":" + id
    ref = elems[3]
    alt = elems[4].split(",")

    # check if this is a multiAllelic
    if remove_multiAllelic and len(alt) > 1:
        return [lineNumber, False, logid + "\tMultiAllelic\t-\n"]

    # check if this variant is an indel
    isIndel = False
    if len(ref) > 1:
        isIndel = True

    if not isIndel:
        for a in alt:
            if len(a) > 1:
                isIndel = True

    filter = elems[6]
    if isIndel:
        # check VQSR first
        if filter.startswith("VQSR") and check_VQSR_Indel:
            if not parseVQSR(filter, isIndel):
                return [lineNumber, False, logid + "\tIndelBelowVQSR\t-\n"]
        elif remove_nonPASS_indel and filter != "PASS":  # if this is an indel, check if it has the PASS label
            return [lineNumber, False, logid + "\tIndelNonPass\t-\n"]
    else:
        # check VQSR first
        if filter.startswith("VQSR") and check_VQSR_SNV:
            if not parseVQSR(filter, isIndel):
                return [lineNumber, False, logid + "\tSNVBelowVQSR\t-\n"]
        elif remove_nonPASS_SNV and filter != "PASS":  # if this is not an indel, check if it has the pass label
            return [lineNumber, False, logid + "\tSNVNonPass\t-\n"]
    # parse info string
    infoelems = elems[7].split(";")
    for elem in infoelems:
        vals = elem.split("=")
        if len(vals) > 1:
            key = vals[0]
            value = vals[1]
            if key == "InbreedingCoeff":
                if value == ".":
                    return [lineNumber, False,
                            logid + "\tIncorrectInbreedingCoeff:" + value + "\t-\n"]
                valuef = float(value)
                if valuef < tresh_IB:
                    return [lineNumber, False,
                            logid + "\tBelowInbreedingCoeff:" + value + "\t-\n"]
        # other filters based on info column can be implemented here...

    elems = line.strip().split("\t")
    # parse genotypes
    format = elems[8].split(":")
    gtcol = -1
    adcol = -1
    dpcol = -1
    gqcol = -1
    for c in range(0, len(format)):
        if format[c] == "GT":
            gtcol = c
        elif format[c] == "AD":
            adcol = c
        elif format[c] == "DP":
            dpcol = c
        elif format[c] == "GQ":
            gqcol = c
    formatOut = ""
    if gtcol > -1:
        formatOut = "GT"
    if dpcol > -1:
        formatOut += ":DP"
    if adcol > -1:
        formatOut += ":AD"
        formatOut += ":AB"
    if gqcol > -1:
        formatOut += ":GQ"

    # only continue if there are genotypes to parse
    if gtcol < 0:
        return [lineNumber, False, logid + "\tNoGTCol\t-\n"]

    sampleinfos = []

    dosagesPreFilter = []
    cctr = 0

    # get current dosages
    # store split columns temporarily
    sampledataelems = [None] * len(elems)
    sampledosages = [-1] * len(elems)
    for i in range(9, len(elems)):
        if includeSample[cctr]:
            # split the sample columns
            sampledata = elems[i].split(":")
            sampledataelems[i] = sampledata
            gt = sampledata[gtcol].split("/")
            if len(gt) < 2:
                # less than two alleles; perhaps the data is phased
                gt = sampledata[gtcol].split("|")
                gtsep = "|"
            if len(gt) < 2 or len(gt) > 2:
                # malformatted genotypes
                # more than two alleles.. skip variant?
                # fhlog.write(headerInfo[1]+"\tPolyPloid?\t-\n"
                dosagesPreFilter.append(-1)
            else:
                if gt[0] == ".":
                    dosagesPreFilter.append(-1)
                else:
                    dosagePrefilter = int(gt[0]) + int(gt[1])
                    dosagesPreFilter.append(dosagePrefilter)
                    sampledosages[i] = dosagePrefilter
        cctr += 1

    # calculate allele frequency etc before filtering
    stats = getStats(dosagesPreFilter, tresh_MAF, tresh_CR, tresh_HWE)
    if not stats[0]:
        return [lineNumber, False,
                logid + "\tFailedPrefilterVarStats\t" + stats[1] + "\n"]

    # check whether variant becomes monomorphic after filtering poor calls
    # iterate samples
    dosagesPostFilter = []
    nrGenotypesReplaced = 0
    averageDepth = 0
    averageDepthCalls = 0
    poorDP = 0
    poorGQ = 0
    poorABHomA = 0
    poorABHomB = 0
    poorABHet = 0
    cctr = 0
    for i in range(9, len(elems)):
        sampledata = sampledataelems[i]
        if sampledata is not None:
            dosage = sampledosages[i]
            ab = -1
            dp = 0
            # if the genotype is not missing
            # check genotype qual
            if dpcol > -1 and dpcol < len(
                    sampledata):  # ugh, VCF really allows all kind of stuff in their genotype format... even if you parse the FORMAT column, there's no telling what is actually provided per sample..
                # determine read depth
                if sampledata[dpcol] == ".":
                    dp = 0
                else:
                    dp = int(sampledata[dpcol])
                    averageDepth += dp
                    averageDepthCalls += 1
                if dp < tresh_DP:
                    dosage = -1
                    poorDP += 1
                    nrGenotypesReplaced += 1
            if gqcol > -1 and dosage > -1 and gqcol < len(sampledata):
                gq = float(sampledata[gqcol])
                if gq < tresh_GQ:
                    if dosage == 0 and gq == 0 and dp >= tresh_DP:
                        # bypass potential error in gVCF merging, see: https://github.com/broadinstitute/gatk/issues/5445
                        # best way would to also parse the PL field, if available, and check whether the homRef and het fields are both 0 (unlikely for a homRef call)
                        pass
                    else:
                        dosage = -1
                        nrGenotypesReplaced += 1
                        poorGQ += 1
            # check allelic balance
            # only need to check this for hets, I guess
            if dosage > -1 and adcol > -1 and adcol < len(sampledata):
                ad = sampledata[adcol].split(",")
                ad1 = float(ad[0])
                ad2 = float(ad[1])
                if ad1 + ad2 == 0:
                    dosage = -1
                    nrGenotypesReplaced += 1
                else:
                    ab = ad2 / (ad1 + ad2)
                    # data is already corrected for DP, so no need to check again.
                    if dosage == 0:
                        if ab > tresh_AB_lower:
                            dosage = -1
                            nrGenotypesReplaced += 1
                            poorABHomA += 1
                    if dosage == 1:
                        if ab < tresh_AB_lower or ab > tresh_AB_upper:
                            dosage = -1
                            nrGenotypesReplaced += 1
                            poorABHet += 1
                    if dosage == 2:
                        if ab < tresh_AB_upper:
                            dosage = -1
                            nrGenotypesReplaced += 1
                            poorABHomB += 1
            dosagesPostFilter.append(dosage)

            sampleinfo = []
            gtout = sampledata[gtcol]
            if dosagePrefilter == -1:
                gtout = "." + gtsep + "."
            sampleinfo.append(gtout)

            # append original information, but note that missing genotypes
            # do not always conform to the FORMAT string
            if dpcol > -1 and dpcol < len(sampledata):
                sampleinfo.append(sampledata[dpcol])
            if adcol > -1 and adcol < len(sampledata):
                sampleinfo.append(sampledata[adcol])
                sampleinfo.append(str(round(ab, 2)))
            if gqcol > -1 and gqcol < len(sampledata):
                sampleinfo.append(sampledata[gqcol])

            sampleinfos.append(":".join(sampleinfo))
        cctr += 1

    # average read depth
    if averageDepthCalls > 0:
        averageDepth /= averageDepthCalls
    else:
        averageDepth = 0

    # variant level quals
    mafPostfilter = getMAF(dosagesPostFilter)
    if mafPostfilter == 0:
        logoutln = (logid + "\tMonomorphicPostFilter\t" + stats[1] +
                    ";NrGenosReplaced:" + str(nrGenotypesReplaced) +
                    ";PoorDP:" + str(poorDP) +
                    ";AvgDP:" + str(round(averageDepth, 2)) + " (" + str(
                    averageDepthCalls) + " calls)" +
                    ";PoorGQ:" + str(poorGQ) +
                    ";PoorABHomA:" + str(poorABHomA) +
                    ";PoorABHomB:" + str(poorABHomB) +
                    ";PoorABHet:" + str(poorABHet) + "\n")
        return [lineNumber, False, logoutln]

    # stats = getStats(dosagesPreFilter, tresh_MAF, tresh_CR, tresh_HWE)

    formatOut = ""
    if gtcol > -1:
        formatOut = "GT"
    if dpcol > -1:
        formatOut += ":DP"
    if adcol > -1:
        formatOut += ":AD"
        formatOut += ":AB"
    if gqcol > -1:
        formatOut += ":GQ"

    elems[8] = formatOut

    if stripINFOcol:
        elems[7] = "."
    outln = "\t".join(elems[0:9])  # construct line header
    # write variant, somehow
    outln += "\t" + "\t".join(sampleinfos) + "\n"
    return [lineNumber, True, logid + "\tPASSQC\t" + stats[1] + "\n", outln]


fh = gzip.open(vcfIn, 'rt')
fho = gzip.open(vcfOut + "-filtered.vcf.gz", 'wt')
fhlog = gzip.open(vcfOut + "-filtered.log.gz", 'wt')
fhlog.write("Id\tReason\tStats\n")

includeSample = []

linectr = 0
lineswritten = 0
for line in fh:
    cctr = 0
    if line.startswith("#CHROM"):
        fho.write(summarizeSettings())
        # header line with samples
        elems = line.strip().split("\t")
        elems[6] = "FILTER"
        elems[7] = "INFO"
        headerout = "\t".join(elems[0:9])
        samplesIncluded = 0
        samplesExcluded = 0
        for i in range(9, len(elems)):
            sample = elems[i]
            if sampleFilterSet is None or sample in sampleFilterSet:
                headerout += "\t" + sample
                includeSample.append(True)
                samplesIncluded += 1
            else:
                includeSample.append(False)
                samplesExcluded += 1
        print("{} samples selected from VCF header, {} excluded".format(
            samplesIncluded, samplesExcluded))
        fho.write(headerout + "\n")
        lineswritten += 1
    elif line.startswith("#"):
        if line.startswith("##FORMAT") and line[13:15] not in ["GT", "DP", "AD", "AB", "GQ"]:
            continue
        if line.startswith("##INFO") and stripINFOcol:
            continue
        fho.write(line)
        lineswritten += 1
    else:
        # this can now be potentially multi-threaded?
        parsed = parseLine(linectr, line,
                           includeSample)  # returns: [lineNumber, False, logid+"\PASSQC\t"+stats[1], outln]
        fhlog.write(parsed[2])
        if parsed[1]:
            fho.write(parsed[3])
            lineswritten += 1
    linectr += 1
    if debug and linectr == stopafterlines:
        print("DEBUG: stopped after " + str(linectr) + " lines")
        break
    if linectr % 10000 == 0:
        print("{} lines parsed, {} written".format(linectr, lineswritten),
              end='\r')
        fho.flush()
        fhlog.flush()
#               if linectr > 100000:
#                       break
print("{} lines parsed, {} written".format(linectr, lineswritten), end='\n')
print("Done. How about that!")
fh.close()
fho.close()
fhlog.close()
