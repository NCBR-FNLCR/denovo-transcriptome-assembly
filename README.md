# Short read based denovo transcriptome assembly pipeline
- based on Snakemake version 5.13.0
- based on Trinity version 2.9.0

# Steps in the pipeline:
1. Creates symbolic links for raw sequencing data in working directory
3. Performs quality control, adapter removal, quality trimming, k-mer based error correction
4. Taxonomic classification of pre-processed read sequences (read contamination)
5. Generate denovo assemblies for each sample separately and all samples combined.
6. Evaluate assembly quality and completeness statistics
7. Cluster assembly contigs to remove redundancy and duplications
8. Search for homologous transcripts against non-redundant nucleotide databases 
9. Gene content and annotation completeness - _still working on it_

# Update the path to working directory (R=) in following files:
- submit_slurm.sh
- pipeline_ctrl.sh 

# Update config.yaml file
- Path to raw data directory
- Path to organism specific lineage BUSCO dataset (https://busco-data.ezlab.org/v4/data/lineages/)
- Taxid ID for organism of interest - homology search
