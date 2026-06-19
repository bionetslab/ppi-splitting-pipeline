# PPI splitting pipeline
Automated leakage-aware splitting of a given PPI dataset into train, validation, and test set.

## Process

Obtain your PPI annotation file `ppis.csv` in the following format: 

```
protein1,protein2
P45985,Q14315
Q86TC9,P35609
...
```

1. The protein sequences `sequences.fasta` are obtained in FASTA format as one-liners (no line breaks):
```
>P32234
MSTILEKISAIESEMARTQKNKATSAHLGLLKAKLAKLRRELISPKGGGGGTGEAGFEVAKTGDARVGFVGFPSVGKSTLLSNLAGVYSEVAAYEFTTLTTVPGCIKYKGAKIQLLDLPGIIEGAKDGKGRGRQVIAVARTCNLIFMVLDCLKPLGHKKLLEHELEGFGIRLNKKPPNIYYKRKDKGGINLNSMVPQSELDTDLVKTILSEYKIHNADITLRYDATSDDLIDVIEGNRIYIPCIYLLNKIDQISIEELDVIYKIPHCVPISAHHHWNFDDLLELMWEYLRLQRIYTKPKGQLPDYNSPVVLHNERTSIEDFCNKLHRSIAKEFKYALVWGSSVKHQPQKVGIEHVLNDEDVVQIVKKV
>P30375
MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFSTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDRNTRNVKAHSQTDRVDLGTLRGYYNQSEDGSHTIQRMYGCDVGSDGRFLRGYQQDAYDGKDYIALNEDLRSWTAADMAAEITKRKWEAAHFAEQLRAYLEGTCVEWLRRHLENGKETLQRTDAPKTHMTHHAVSDHEAILRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPEPLTLRWEPSSQPTIPIVGIIAGLVLFGAVIAGAVVAAVRWRRKSSDRKGGSYSQAASSDSAQGSDVSLTACKV
>P16209
MAVMPPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDEETRSAKAHSQTDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHAAEQRRAYLEGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSLTACKV
>P16211
MAIMAPRTLLLLLSGALALTQTWAGSHSMRYFSTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRTPWMEQEGPEYWDRETRSVKAHAQTNRVDLGTLRGYYNQSDGGSHTIQRMFGCDVGPDGRFLRGYEQHAYDGKDYIALNEDLRSWTAADMAAQITQRKWEAAGAAEQDRAYLEGLCVESLRRYLENGKETLQRTDAPKTHMTHHPVSDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGKEQRYTCHVQHEGLPEPLTLRWELSSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRRNSDRKGGSYSQAASNDSAQGSDVSLTACKV
```

2. Sequence lengths are obtained for each protein for a length-normalized bitscore


3. All-against-all similarities are calculated with BLAST:

```{bash}
makeblastdb -dbtype prot -in sequences.fasta
blastp -query sequences.fasta -db mydb -outfmt "6 qseqid sseqid evalue bitscore"  -max_hsps 1 -out all_vs_all.tsv
```

4. Make a similarity network and save it as METIS file for `KaHIP`. Edge weights can be specified as bitscore or length-normalized bitscore.

5. Install `KaHIP locally`

6. Run `KaHIP's kaffpa` command on the METIS file

```{bash}
kaffpa similarity.graph --seed=1234 --output_filename="partitioned_proteome.txt" --k=3 --preconfiguration=strong
```

7. Sort the PPIs into their partition blocks. PPIs with proteins belonging to different blocks are discarded. The largest resulting block is used as the training set, the second largest as the validation set, and the smallest as the test set.

8. Three sequence files are written for the proteins belonging to the training, validation, and test set, respectively. 

9. The sequence files are used to run CD-HIT to reduce redundancy between the blocks:

```{bash}
cd-hit-2d -i Intra_0.fasta -i2 Intra_1.fasta -o sim_intra0_intra_1.out -c 0.4 -n 2
cd-hit-2d -i Intra_0.fasta -i2 Intra_2.fasta -o sim_intra0_intra_2.out -c 0.4 -n 2
cd-hit-2d -i Intra_1.fasta -i2 Intra_2.fasta -o sim_intra1_intra_2.out -c 0.4 -n 2
```

10. The redundant sequences are removed from the training, validation, and test set, respectively.

11. Negatives are sampled randomly for each block separately such that, e.g., all proteins occurring in training negative samples occur in the positive training samples, and so on. The negatives are sampled in a way that, in expectation, preserves the individual node degrees of positive samples, such that a protein has approximately the same number of positive and negative annotations.

