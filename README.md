# ProteinBERT for Huggingface JAX/Flax community week
## 1. PreTrain ProteinBERT from scratch
ProteinBERT is a universal deep-learning model of protein sequence and function based on the BERT architecture. The goal of this project is to pretrain the ProteinBERT model in JAX/Flax for downstream finetuning tasks like predicting protein structure, post translational modifications and/or biophysical attributes. The model itself can be seen as an extension to the Protein Transformer models trained by the Rostlab.

## 2. Language
The model will be trained on protein sequences and gene ontology (GO) annotations.

## 3. Model
The classic Transformer/BERT architecture, with some additions.

## 4. Datasets
ProteinBERT is pretrained on a dataset derived from UniRef90 which consists of ~106M protein sequences. Protein Gene Ontology (GO) annotations are used as additional inputs (to help the model infer about the function of the input protein and update its internal representations and outputs accordingly).

Possible links to publicly available datasets include:
- [GO annotations from CAFA](https://www.biofunctionprediction.org/cafa-targets/cafa4ontologies.zip)
- [UniRef90 XML](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.xml.gz)
- [UniRef90 FASTA](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz)

## 5. Training scripts
Data preprocessing scripts are available on the [Official Github Repo](https://github.com/nadavbra/protein_bert). The training script can be adapted from [run_mlm_flax.py](https://github.com/huggingface/transformers/blob/master/examples/flax/language-modeling/run_mlm_flax.py).

## 6. Challenges
- The data of protein sequences and GO annotations require ~1 TB of scratch disk space.
- The original paper states a pretraining time of ~28 days on a single GPU (Nvidia Quadro RTX 5000).
- Besides NLP skills there may be some bioinformatics skills required

## 7. Desired project outcome
A model that can be further finetuned to predict protein structure, post translational modifications, and/or biophysical attributes.

## 8. Reads
The following links can be useful to better understand the project and what has previously been done.
- [ProteinBERT article](https://www.biorxiv.org/content/10.1101/2021.05.24.445464v1)
- [Official Github Repo](https://github.com/nadavbra/protein_bert)
- [Work in progress porting ProteinBERT to JAX by Saurav Maheshkar](https://github.com/SauravMaheshkar/ProteinBERT)
