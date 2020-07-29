max_hamming_dist=1
indrops_pipeline_outdir=Hashtags_lead/
outdir=Hashtags_lead/output/hash${max_hamming_dist}

mkdir -p $outdir
wd=$PWD
cd ${indrops_pipeline_outdir}

for f in *; do
    echo $f
    cat $f/filtered_parts/*sorted.fastq.gz | gunzip -c | python $wd/get_cell_hashing_barcodes_flex.py $wd/expected_barcodes.txt $wd/$outdir/$f.counts.csv $wd/$outdir/$f.reads.csv $wd/$outdir/$f.pickle ${max_hamming_dist}
    echo ""
done
