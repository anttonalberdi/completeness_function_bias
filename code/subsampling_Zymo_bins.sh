# Raphael Eisenhofer 2023

Zymo reads (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9223867/)

conda activate kingfisher

kingfisher get -r SRR12324251 -m aws-http --download-threads 16

seqtk sample -s 1337 SRR12324251_1.fastq 1000000 > SRR12324251-1Mss_1.fastq
seqtk sample -s 1337 SRR12324251_2.fastq 1000000 > SRR12324251-1Mss_2.fastq

seqtk sample -s 1337 SRR12324251_1.fastq 2000000 > SRR12324251-2Mss_1.fastq
seqtk sample -s 1337 SRR12324251_2.fastq 2000000 > SRR12324251-2Mss_2.fastq

seqtk sample -s 1337 SRR12324251_1.fastq 5000000 > SRR12324251-5Mss_1.fastq
seqtk sample -s 1337 SRR12324251_2.fastq 5000000 > SRR12324251-5Mss_2.fastq

pigz -p 16 *.fastq

# Ran the EHI pipeline: https://github.com/earthhologenome/EHI_bioinformatics
# 1_Preprocess_QC.snakefile
# 2_Individual_Assembly_Binning.snakefile


for bin in *.fa.gz; 
    do for seed in `seq 1 10`;
       do for rate in `seq 0.6 0.1 0.9`;
          do reformat.sh \
                in=$bin \
                out=${bin/.fa.gz/_SEED"$seed"_RATE"$rate".fa.gz} \
                sampleseed=$seed \
                samplerate=$rate;
            done;
    done;
done
