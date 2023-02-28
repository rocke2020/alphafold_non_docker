import argparse, os
from datetime import datetime


DATE_TIME = "%Y_%m_%d %H:%M:%S"


class ArgparseUtil(object):
    """
    参数解析工具类
    """
    def __init__(self):
        """ Basic args """
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("--seed", default=2, type=int)
        self.parser.add_argument('--gpu_device_id', default=2, type=int, help='the GPU NO.')

    def classifier_cnn_attention(self):
        """ task args """
        self.parser.add_argument("--model_type", type=str, default='esm', help="prot or esm as pretrained embeddings")
        self.parser.add_argument("--task_name", type=str, default='classifier', help="classifier")
        self.parser.add_argument("--enable_train", type=int, default=1, help="0 is false, 1 is true")
        self.parser.add_argument("--enable_final_train", type=int, default=0,
            help="split data only 2 parts, train and val; 0 is false, 1 is true")
        self.parser.add_argument("--pre_train", type=int, default=0, help="test auto lr finder, etc. 0 false, 1 true")
        self.parser.add_argument("--embedding_type", type=str, default='', help="vocab, esm, t5")
        self.parser.add_argument("--use_fp16_in_pretrained_embedding", type=int, default=1, help="bool")
        self.parser.add_argument("--read_dataset_cache", type=int, default=1, help="0 is false; 1 true")
        self.parser.add_argument("--enable_over_sample", type=int, default=1, help="0 is false; 1 true")
        self.parser.add_argument("--cache_path", type=str, default='CAMP/cache', help="")
        self.parser.add_argument("--lightning_logs_root_path", type=str, default='.', help="")
        self.parser.add_argument('--model_path', type=str, help='the path of trained model path for predict',
            default="CAMP/cache/lightning_logs/version_24/checkpoints/epoch=09-test_auprc=0.8358-test_auc=0.9497-test_acc=0.9067-test_loss=0.30192-train_loss=0.0787.ckpt",  # v2.3
        )
        self.parser.add_argument('--protein_fasta_path', dest='protein_fasta_path',
            default="CAMP/data/protein.fasta", type=str, help='the path of protein.fasta')
        self.parser.add_argument('--peptide_fasta_path', dest='peptide_fasta_path',
            default="CAMP/data/peptide.fasta", type=str, help='the path of peptide.fasta')
        self.parser.add_argument('--pretrain_model_path', dest='pretrain_model_path',
            default="/mnt/sda/models/esm/checkpoints/esm2_t33_650M_UR50D.pt", type=str, help='')

        self.parser.add_argument("--enable_cluster", type=int, default=1, help="0 is false, 1 is true")
        self.parser.add_argument("--selected_cluster_nums", type=list,
            default=[3], help="-1 means run train on all clusters; 0~n, measn the cluter num")

        self.parser.add_argument('--data_split_seed', default=0, type=int, help='the seed to split full data')
        self.parser.add_argument('--data_version', default='2.5', type=str, help='')

        # training hyperparameter
        self.parser.add_argument('--batch_size', default=32, type=int)
        self.parser.add_argument('--val_batch_size', default=32, type=int)
        self.parser.add_argument('--n_fold', dest='n_fold',
            default=5, type=str, help='the number of cross validation fold')
        self.parser.add_argument('--epochs', default=60, type=int, help='the number of epochs')
        self.parser.add_argument('--learning_rate', default=0.00009, type=str, help='')
        self.parser.add_argument('--momentum', default=0.9, type=str, help='the number of momentum')
        self.parser.add_argument('--wd', default=0.0001, type=str, help='the number of weight decay')
        
        self.parser.add_argument("--dropout", type=float, default=0.2, help="")
        self.parser.add_argument("--warmup_steps", type=float, default=0.0,
            help="if warmup_steps < 1, it is used as int(warmup_ratio * max_interation_num)")

        self.parser.add_argument("--enable_extra_predict", type=int, default=1, help="0 is false, 1 is true")
        self.parser.add_argument("--predict_input_file_type", type=str,
            default='interpep',
        )

        args = self.parser.parse_args()
        return args


def save_args(args, output_dir='.', with_time_at_filename=False):
    os.makedirs(output_dir, exist_ok=True)

    t0 = datetime.now().strftime(DATE_TIME)
    if with_time_at_filename:
        out_file = os.path.join(output_dir, f"args-{t0}.txt")
    else:
        out_file = os.path.join(output_dir, f"args.txt")
    with open(out_file, "w", encoding='utf-8') as f:
        f.write(f'{t0}\n')
        for arg, value in vars(args).items():
            f.write(f"{arg}: {value}\n")


def log_args(args, logger):
    for arg, value in vars(args).items():
        logger.info(f"{arg}: {value}")
