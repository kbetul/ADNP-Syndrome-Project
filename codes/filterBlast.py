iFile = open('../data/ADNP_blast_result1000.fas', 'r')
oFile = open('../data/ADNP_filtered_isoform1000.fas', 'w')

skip = True

for line in iFile:
    if line.startswith('>'):
        if 'isoform' in line:
            skip = True
        else:
            skip = False
            oFile.write(line)
    elif not skip:
        oFile.write(line)
oFile.close()
iFile.close()    