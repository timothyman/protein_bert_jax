Processing pretraining data for ProteinBERT
=============

We warn that this is a long process and it also requires a lot of storage (>1TB).

Step 1: Create the UniRef dataset
------------

ProteinBERT can be pretrained on a dataset derived from either UniRef50 (recommended) or UniRef90. Follow these steps to produce this dataset:

1. First, choose a working directory with sufficient (>1TB) free storage.

.. code-block:: sh
    
    cd /some/workdir

2. Download the metadata of GO from CAFA and extract it.

.. code-block:: sh

    wget https://www.biofunctionprediction.org/cafa-targets/cafa4ontologies.zip
    mkdir cafa4ontologies
    unzip cafa4ontologies.zip -d cafa4ontologies/
    
3. Download UniRef50 as both XML, substitute 50 with 90 for UniRef90.

.. code-block:: sh

    wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.xml.gz

4. Use the *create_uniref_h5_dataset_onestep* script to create the final dataset (in the H5 format). An intermediate step is taken to extract the GO annotations associated with UniRef's records into an SQLite database (and a CSV file with the metadata of these GO annotations). Since this is a long process (which can take up to a few days), it is recommended to run this in the background (e.g. using *nohup*).

.. code-block:: sh
    
    nohup ./create_uniref_h5_dataset_onestep --uniref-xml-gz-file=./uniref50.xml.gz --go-annotations-meta-file=./cafa4ontologies/go.txt --intermediate-sqlite-db-file=./uniref50_proteins_and_annotations.db --intermediate-go-annotations-meta-csv-file=./uniref50_go_annotations.csv --output-h5-dataset-file=./uniref500_dataset.h5 --min-records-to-keep-annotation=100 >& ./log_create_uniref50_h5_dataset.txt &
    
    
5. If you are planning to evaluate your model on certain downstream benchmarks, it is recommended that any UniRef record similar to a test-set protein in these benchmark will be considered part of the pretraining's test set. You can use provide them to *set_h5_testset*  (step 6) through the flag ``--uniprot-ids-file=./uniref_50_seqs_matching_test_set_seqs.txt``, where the provided text file contains the UniProt IDs of the relevant records, one per line (e.g. *A0A009EXK6_ACIBA*). You can use the *create_test_ids* script to crawl through your UniRef50 xml database again and through your downstream benchmarks files (assumption is that test datasets are ``*.test.csv`` formats) as saved in the test folder. 

.. code-block:: sh
    
    nohup ./create_test_ids --uniref-xml-gz-file=./uniref50.xml.gz --test-file-folder=../test_data --output-file=./uniref_50_seqs_matching_test_set_seqs.txt >& ./log_create_test_ids_.txt &
    
6. Finally, use ProteinBERT's *set_h5_testset* script to designate which of the dataset records will be considered part of the test set (so that their GO annotations are not used during pretraining). 

.. code-block:: sh

    set_h5_testset --h5-dataset-file=./uniref50_dataset.h5 --uniprot-ids-file=./uniref_50_seqs_matching_test_set_seqs.txt

