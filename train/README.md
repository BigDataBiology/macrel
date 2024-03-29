# Code for training Macrel models

## DATA

To run this code, you must download the data from the original sources:

1. _AmPEP: Sequence-based prediction of antimicrobial peptides using
   distribution patterns of amino acid properties and random forest_ by Pratiti
   Bhadra, Jielu Yan, Jinyan Li, Simon Fong &amp; Shirley W. I. Siu Scientific
   Reports volume 8, Article number: 1697 (2018)
   https://www.nature.com/articles/s41598-018-19752-w#Sec9
2. _iAMP-2L: a two-level multi-label classifier for identifying antimicrobial
   peptides and their functional types_ Xiao X, Wang P, Lin WZ, Jia JH, Chou
   KC. Anal Biochem. 2013, 436(2):168-77. doi: 10.1016/j.ab.2013.01.019.
3. _A Web Server and Mobile App for Computing Hemolytic Potency of Peptides_
   Kumardeep Chaudhary, Ritesh Kumar, Sandeep Singh, Abhishek Tuknait, Ankur
   Gautam, Deepika Mathur, Priya Anand, Grish C. Varshney &amp; Gajendra P. S.
   Raghava Scientific Reports volume 6, Article number: 22843 (2016)
   https://www.nature.com/articles/srep22843

You can adapt the code below:

```bash

mkdir -p data
cd data/

wget https://github.com/ShirleyWISiu/AmPEP/raw/master/M_model_train_AMP_sequence.zip
wget https://github.com/ShirleyWISiu/AmPEP/raw/master/M_model_train_nonAMP_sequence.zip
unzip M_model_train_AMP_sequence.zip
unzip M_model_train_nonAMP_sequence.zip

wget -O hemo.training.pos.faa 'https://webs.iiitd.edu.in/raghava/hemopi/data/HemoPI_1_dataset/main/pos.fa'
wget -O hemo.training.neg.faa 'https://webs.iiitd.edu.in/raghava/hemopi/data/HemoPI_1_dataset/main/neg.fa'

wget -O hemo.validation.pos.faa https://webs.iiitd.edu.in/raghava/hemopi/data/HemoPI_1_dataset/validation/pos.fa
wget -O hemo.validation.neg.faa https://webs.iiitd.edu.in/raghava/hemopi/data/HemoPI_1_dataset/validation/neg.fa


wget http://www.jci-bioinfo.cn/iAMP/Supp-S1.pdf
wget http://www.jci-bioinfo.cn/iAMP/Supp-S2.pdf
```

Then you must **extract the dataset from `Supp-S1.pdf` and `Supp-S2.pdf`** and
save it as four FASTA files: `Supp-S[12]_AMP.faa` and `Supp-S[12]_NAMP.faa`
(`pdftotext` can partially automate this, but it still requires a bit of manual
labour).

You should have 8 files in the end (inside the `data/` subdirectory):

- `M_model_train_AMP_sequence.fasta`
- `M_model_train_nonAMP_sequence.fasta`
- `Supp-S2_AMP.faa`
- `Supp-S2_NAMP.faa`
- `hemo.training.neg.faa`
- `hemo.training.pos.faa`
- `hemo.validation.neg.faa`
- `hemo.validation.pos.faa`

## Training/testing models

Training the models is a two-step process:

1. Generate the feature tables
2. Build the models & testing

Note that `macrel` needs to be installed as the Python scripts import it

```bash
python build-AMP-table.py
python build-Hemo-table.py

python train-models.py
```

The `train-models.py` script will also perform testing. In the case of the AMP
tests, it does it in a careful way as it needs to avoid having the same
sequences in the training and testing sets.

