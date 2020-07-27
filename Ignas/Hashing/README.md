# Single-cell Hashing
Single-cell hashing methods are a group of molecular tools designed for sample multiplexing in droplet-based single cell transcriptome barcoding platforms. They are based on genetic labeling of cells with sample-specific barcodes before pooling and single-cell isolation. These genetic barcodes possess a poly-A tail rendering their capture during transcriptome barcoding reaction. Hashing approaches address the constrains of parallel processing of individual samples and expand the limitation of Poissonian statistics based single-cell encapsulation events. Both of these factors affect the sample processing throughput and increase the consumption of reagents and time.
[Scheme](/misc/scheme1.png)
# Extracting the hashtag barcodes
When hashing is used with InDrops, molecules that house sample index, cell barcode and hashtag (HTO) barcode are created. 
1.	Run Standard InDrops pipeline on the transcriptome data (If you don’t want the Transcriptome data and you only did MiSeq on the hashtags, ignore this step).
2.	First tree steps of InDrops pipeline can be used to generate FASTQ files that house the sample index, cell barcode and HTO sequence.
    - This follows the same procedure as InDrops pipeline run. You set up a .yaml file (check [my_analysis_hashtags.yaml](/my_analysis_hashtags.yaml)). It is identical to the original InDrops.yaml file, except I have included a Trimmomatic argument CROP and complexity filtering threshold was reduced. The HTO barcode sequences are short, therefore it is better to trim the reads to a shorter length to avoid losing data due to adapter sequences. You do not need to change the reference as these reads are not used for aligning.
    - **NOTE**: If MiSeq is used, it is better to include additional trimmomatic arguments as reads tend to be lower quality:
      - trimmomatic_arguments:
        - LEADING: "0"
        - SLIDINGWINDOW: "0:0"
        - MINLEN: "0"
        - CROP: "25"
        - argument_order: ['LEADING', 'SLIDINGWINDOW', 'MINLEN', 'CROP']
    - Submit InDrops job only for the first 3 steps (check [submit_indrops_hashtags.sh](/submit_indrops_hashtags.sh)) and this creates an output folder.
3.	This directory is then used as input for a custom python code (**wrote by Sam Wolock**) to generate cells x HTO count matrices (in .csv format).
    - Add the hashtag barcode sequences to the [expected_barcodes.txt](/expected_barcodes.txt) file.
    - Edit the [run_cell_hashing_barcodes_clean.sh](/run_cell_hashing_barcodes_clean.sh) file. You need to change the **indrops_pipeline_outdir** to the directory of the 1st step output. Also adjust the **hamming distance** as needed.
    - All the files (expected_barcodes.txt, run_cell_hashing_barcodes_clean.sh, get_cell_hashing_barcodes.py and get_cell_hashing_barcodes_flex.py) have to be in the same folder. 
    - Start an interactive job on the O2 and bash the **run_cell_hashing_barcodes_clean.sh** file. 
    - This creates an output folder that houses the cells x HTO counts, cells x HTO reads .csv and a pickle files for each of the library in the InDrops pipeline output folder
4.	Use the output with the python notebook to analyze the data.
  

