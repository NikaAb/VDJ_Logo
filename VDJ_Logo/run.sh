python3 biologo.py --logotype -mv --alignments -t Dataset/logo22_07/PIOL.txt Dataset/logo22_07/PCNSL.txt Dataset/logo22_07/DLBCL.txt Dataset/logo22_07/Nonsubsets_CLL.txt Dataset/logo22_07/subsets_CLL.txt Dataset/logo22_07/MCL.txt --germline -t Dataset/logo22_07/IGVH4_34.txt --labels PIOL PCNSL DLBCL  CLLNonSubsets CLLSubsets MCL --output -p MV_22_7 --settings -seqtype p -start 1 -length -1 -blocksize 200
python3 biologo.py --logotype -p --alignments -t  Dataset/logo22_07/PIOL.txt --germline -t Dataset/logo22_07/IGVH4_34.txt  --output -p SL --settings -seqtype p -start 1 -length -1 -blocksize 200
python3 biologo.py --logotype -s --alignments -f  Dataset/subfamilylogo/family.aln --germline -f Dataset/subfamilylogo/subfamily.aln  --output -p FL --settings -seqtype p -start 1 -length -1 -blocksize 200
python3 biologo.py --logotype -mg --input mutagenesis/mutagenesis.txt  --output -p MG --settings -seqtype p -start 1 -length -1 -blocksize 50 -posth -0.5 -negth -6