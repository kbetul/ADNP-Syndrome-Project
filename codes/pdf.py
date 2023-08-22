import fitz
from pprint import pprint
from findConservation import *
from progress.bar import Bar

#Code taken from https://github.com/pymupdf/PyMuPDF/discussions/1532
def changeColorText(text : str, doc : fitz, txtColor : tuple) -> None:
    # text = "NP_001269460.1_activity-dependent_neuroprotector_homeobox_protein_Homo_sapiens"
    # doc = fitz.open("./data/ADNP_tree1000.pdf")
    page = doc[0]
    rl = page.search_for(text)
    assert len(rl) == 1
    clip = rl[0]
    # extract text info now - before the redacting removes it.
    blocks = page.get_text("dict", clip=clip)["blocks"]
    span = blocks[0]["lines"][0]["spans"][0]
    assert span["text"] == text

    # remove text
    page.add_redact_annot(clip)
    page.apply_redactions()

    # re-insert same text - different color
    font = fitz.Font("helvetica-bold")  # this must be known somehow - or simply try some font else
    tw = fitz.TextWriter(page.rect, color=txtColor)
    # text insertion must use the original insertion poin and font size.
    # if not original font, then some fontsize adjustments will probably be required:
    # check bbox.width against textlength computed with new font
    tw.append(span["origin"], text, font=font, fontsize=span["size"])
    tw.write_text(page)
    # doc.ez_save("./outputs/colored_tree1000.pdf")

def colorPDFForVUS(aminoAcid : str, numOfHits : int, numOfVUS : int, msa : dict, index : int, type : str):
    blosum = createBlosum62Dict()
    colors = {
            'green' : (0,255,0),
            'red' : (255,0,0),
            'blue' : (0,0,255)
        }
    color = (0,0,0)
    doc = fitz.open(f"./data/ADNP_tree{numOfHits}.pdf")

    bar2 = Bar(f'Processing Species for VUS{numOfVUS} for {type}', max=len(msa))

    for specie, seq in msa.items():
        newSpecie = specie.replace('[','').replace(']','')

        if 'PREDICTED:' not in specie and 'LOW QUALITY PROTEIN:' not in specie:
            newSpecie = newSpecie.replace(' ', '_')
        
        color = (0,0,0)
        if seq[index] != '-':
            if seq[index] == aminoAcid:
                color = colors['green']
            elif blosum[aminoAcid][seq[index]] > 0:
                color = colors['blue']
            elif blosum[aminoAcid][seq[index]] < 0:
                color = colors['red']
        
        changeColorText(newSpecie, doc, color)
        bar2.next()
    bar2.finish()
    doc.ez_save(f"./colored_trees/colored_tree{numOfHits}_VUS{numOfVUS}_{type}.pdf")

vusDict = createVUSDict()
msa = msaReader(f"./data/ADNP_filtered_isoform1000_msa_muscle.fas")
humanSeq = list(msa.values())[0]

for index in range(len(vusDict)):
    origIndex = list(vusDict.keys())[index]

    msaIndex = findMSAIndex(list(msa.values())[0], origIndex)
    colorPDFForVUS(vusDict[origIndex]['original'], 1000, index + 1, msa, msaIndex - 1, 'original')
    colorPDFForVUS(vusDict[origIndex]['subsitution'], 1000, index + 1, msa, msaIndex - 1, 'substitution')