input_dir=/mnt/sdc/af_input
# input file: 1a1r_D_B.fasta 1a1r_A.fasta
cd alphafold-2.3.1
nohup ./run_alphafold.sh \
    -m multimer \
    -f $input_dir/1a1r_D_B.fasta \
    -d /home/qcdong/af_data \
    -o /mnt/sdc/af_out \
    -t 2022-01-01 \
    -a 1
    > ../zlog/run.log 2>&1 &
