def msaReader(file_name : str) -> dict:
    inputFile = open(file_name, 'r')
    msaDict = {}
    seqName = ""

    for line in inputFile:
        if line.startswith('>'):
            seqName = line[1:].strip()
            msaDict[seqName] = ""
        else:
            msaDict[seqName] += line.strip()
    
    return msaDict

def createBlosum62Dict() -> dict:
    # source = https://gist.github.com/davidball/69e451a8b7e693337988fdb1997bd61f
    blosum62txt = """A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1"""

    rows = [line.split(' ') for line in blosum62txt.split('\n')]
    rows = [[score for score in line if score != ''] for line in rows]
    aminoAcidList = [line[0] for line in rows]
    
    blosumDict = {}
    aminoAcid = ''

    for line in rows:
        counter = 0
        for score in line:
            if score in aminoAcidList:
                blosumDict[score] = {}
                aminoAcid = score
            else:
                blosumDict[aminoAcid][aminoAcidList[counter]] = int(score)
                counter += 1

    return blosumDict

def findMSAIndex(humanSeq : str, index : int) -> int:
    msa_index = 0
    counter = index

    for char in humanSeq:
        if char != '-':
            counter -= 1
        
        msa_index += 1

        if counter == 0:
            break
    return msa_index

def findConservationRatio(msa : dict, index : int, originalAminoAcid : str, substituedAminoAcid : str):
    identityRatioSub = 0
    identityRatioOrig = 0
    positiveRatioSub = 0
    positiveRatioOrig = 0
    zeroRatioSub = 0
    zeroRatioOrig = 0
    negativeRatioSub = 0
    negativeRatioOrig = 0
    gapRatioSub = 0
    gapRatioOrig = 0

    blosumDict = createBlosum62Dict()

    humanSeq = list(msa.values())[0]

    msa_index = findMSAIndex(humanSeq, index)    

    for seq in msa.values():
        if seq[msa_index] == substituedAminoAcid:
            identityRatioSub += 1
        if seq[msa_index] == '-':
                gapRatioSub += 1
        elif blosumDict[seq[msa_index]][substituedAminoAcid] > 0:
            positiveRatioSub += 1
        elif blosumDict[seq[msa_index]][substituedAminoAcid] == 0:
            zeroRatioSub += 1
        else:
            negativeRatioSub += 1
        
        if seq[msa_index] == originalAminoAcid:
            identityRatioOrig += 1
        if seq[msa_index] == '-':
            gapRatioOrig += 1
        elif blosumDict[seq[msa_index]][originalAminoAcid] > 0:
            positiveRatioOrig += 1
        elif blosumDict[seq[msa_index]][originalAminoAcid] == 0:
            zeroRatioOrig += 1
        else:
            negativeRatioOrig += 1
    
    resultDict1, resultDict2 = {}, {}

    identityRatioOrig = float(identityRatioOrig)/len(msa)
    resultDict1['identity'] = identityRatioOrig
    positiveRatioOrig = float(positiveRatioOrig)/len(msa)
    resultDict1['positive'] = positiveRatioOrig
    zeroRatioOrig = float(zeroRatioOrig)/len(msa)
    resultDict1['zero'] = zeroRatioOrig
    negativeRatioOrig = float(negativeRatioOrig)/len(msa)
    resultDict1['negative'] = negativeRatioOrig
    gapRatioOrig = float(gapRatioOrig)/len(msa)
    resultDict1['gap'] = gapRatioOrig

    identityRatioSub = float(identityRatioSub)/len(msa)
    resultDict2['identity'] = identityRatioSub
    positiveRatioSub = float(positiveRatioSub)/len(msa)
    resultDict2['positive'] = positiveRatioSub
    zeroRatioSub = float(zeroRatioSub)/len(msa)
    resultDict2['zero'] = zeroRatioSub
    negativeRatioSub = float(negativeRatioSub)/len(msa)
    resultDict2['negative'] = negativeRatioSub
    gapRatioSub = float(gapRatioSub)/len(msa)
    resultDict2['gap'] = gapRatioSub

    return resultDict1, resultDict2

def createVUSDict() -> dict:
    fileName = './data/vus.txt'
    inputFile = open(fileName, 'r')

    vusDict = {}

    for line in inputFile:
        line = line.strip().split(' ')
        vusDict[int(line[0])] = {'original' : line[1], 'subsitution' : line[2]}
    
    return vusDict

def main():
    file_name = 'ADNP_filtered_isoform1000_msa_muscle.fas'
    msa_dict = msaReader(file_name)
    vusDict = createVUSDict()
    blosum = createBlosum62Dict()

    outFile = open('VUSMetrics.txt', 'w')
    counter = 0

    for index, vus in vusDict.items():
        counter += 1
        original, substitution = findConservationRatio(msa_dict, index-1, vus['original'], vus['subsitution'])

        outFile.write(f"VUS#{counter}\n")
        outFile.write(f"\nOriginal Amino Acid {vus['original']} substitued to {vus['subsitution']} which has BLOSUM62 score of {blosum[vus['original']][vus['subsitution']]}\n")
        outFile.write("\nOriginal Values\n")

        for ratio, score in original.items():
            outFile.write(str(ratio) + " " + str(score) + "\n")

        outFile.write("\nSubstitution Values\n")

        for ratio, score in substitution.items():
            outFile.write(str(ratio) + " " + str(score) + "\n")

        outFile.write("\n")       

    outFile.close()

if __name__=='__main__':
    main()