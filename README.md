# VDJ logo for BCR/TCR intraclonal diversity analysis 
There are currently two separate folders that are not connected in any form. The idea is to create a pipeline using the programs of the "Alignment_clone_sequences" folder and then the "VDJ_Logo" folder, to build the logo from an MSA file containing BCR sequences obtained from clonally related cells as well as the hypothetical or the real naive BCR sequence, as the unmuted version of BCRs within a clone, used as the reference. The different steps of the pipeline would be :

1- Run the "alignment_seq.py" script from the "Alignment_clone_sequences" folder.
2- Using the output generated by the previous step, run the "biologo.py" of the "VDJ_Logo" folder.

Here are some details about each script.


- alignment_seq.py

```
$python alignment_seq.py -a AIRR.tsv -n CloneName
```

Where AIRR.tsv is the annotated BCR/TCR sequences in AIRR format, and CloneName is the output name that will be used to generate 6 files, 3 of which will be used for the next step. These files are the following.

CloneName_germline.fasta : the naive sequence used as the refrence 


CloneName_selected_seq_uniq.aln.fa : the alignment of BCR/TCR sequences


CloneName_region.txt : contning the beginging and end posion of different regions (V, J, cdr1, cdr2, cdr3, D) of each sequence.




The other three files can be deleted at the end of this step.

To run the test, please use the folowing command :

```
$python alignment_seq.py -a Input/sc1_BORJ_clone2_100_biologo.tsv -n test
```



- biologo.py

```
$python biologo.py --logotype [-p|-m|-s|-mv|-mg] --alignments [-f|-a|-t] filename --germline [-f|-a|-t] filename --output [-e|-s|-p] filename --settings -seqtype p -start 1 -length -1 -blocksize 200
```

--alignments takes the CloneName_selected_seq_uniq.aln.fa file and --germline takes CloneName_germline.fasta, for more details please check the biologo.py file. 

If you want to generate the logo for a segment of the alignment (for instance the CDR3 or the V gene, etc) you can use the position information in the CloneName_region.txt file and run the command with the corresponding -start and -length option.

To run the test, please use the folowing command :
```
$python biologo.py  --logotype -s --alignments -t Dataset/test2.txt   --germline -t  Dataset/test1.txt --output -p FL --settings -seqtype p -start 1 -length -1 -blocksize 200
```
