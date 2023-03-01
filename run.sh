[ ! -d zlog ] && mkdir -p zlog
input_dir=/mnt/sdc/af_input
# input file: 1a1r_B_D.fasta 1a1r_A.fasta
input_file=6i51_H_I_pdb_fasta_pos_human.fasta
cd alphafold-2.3.1
nohup ./run_alphafold.sh \
    -m multimer \
    -f $input_dir/$input_file \
    -d /home/qcdong/af_data \
    -o /mnt/sdc/af_out \
    -t 2022-01-01 \
    -a 1 \
    > ../zlog/$input_file.log 2>&1 &
