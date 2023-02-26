input_dir=/mnt/sdc/af_input
cd alphafold-2.3.1
nohup ./run_alphafold.sh \
    -d /home/qcdong/af_data \
    -o /mnt/sdc/af_out \
    -f $input_dir/1a1r_A.fasta \
    -t 2022-01-01 \
    > run.log 2>&1 &
