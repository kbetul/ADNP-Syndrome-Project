# ENS210 Project Group 3 Proposal
 
Upon our research, we came upon a rare Mendelian disease called <a href = "https://rarediseases.info.nih.gov/diseases/12931/adnp-syndrome">ADNP Syndrome</a>. This syndrome is caused by a gene truncation event on the ADNP gene which is one of the causes of the behavioral disorder Autism Spectrum Disorder. The <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6487024/"> study done by Bend et al.</a> showcases the ability to correctly predict undiagnosed patients with ADNP Syndrome according to the model created by the identification of DNA methylation episignatures that indicate ADNP mutation. The study also claims the existence of variants of unknown significance (VUS) that we found examples of through <a href = "https://databases.lovd.nl/shared/variants/ADNP/unique#object_id=VariantOnTranscriptUnique%2CVariantOnGenome&id=ADNP&order=VariantOnTranscript%2FDNA%2CASC&search_transcriptid=00024223&search_VariantOnGenome/ClinicalClassification=VUS&page_size=100&page=1">LOVD (Leiden Open Variation Database)</a> the exact sequence for and the specific substitution that happened. To determine the effect of these variants, we are going to p-BLAST the sequence of the healthy ADNP gene to find homologous proteins in other species. Then we will perform a multiple-sequence alignment in MEGA11. The resulting multiple protein sequences will act as a basis for us to determine whether a VUS is a pathogenic substitution or an allowed substitution. We will reach these conclusions based on the degree of preservation of the amino acid or the substitution of chemically compatible amino acids in that position among the subject protein sequences. MEGA11 allows us to see chemical similarities by annotating amino acids with the same color. For instance, if an overwhelming majority of homologous proteins have a certain or chemically similar amino acid conserved throughout with little to no instances that allowed the substitution observed in the VUS then we conclude that this is a pathogenic substitution.