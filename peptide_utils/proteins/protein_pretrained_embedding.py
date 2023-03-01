import torch
import esm
import logging, os, sys
from functools import cache
from torch.nn.utils.rnn import pad_sequence
from transformers import T5Tokenizer, T5EncoderModel, BertModel, BertTokenizer
sys.path.append(os.path.abspath('.'))
from utils.train_util import get_device
import re


def get_logger(name=None, log_file=None, log_level=logging.DEBUG):
    """ default log level DEBUG """
    logger = logging.getLogger(name)
    logging.basicConfig(format=fmt, datefmt=datefmt)
    if log_file is not None:
        log_file_folder = os.path.split(log_file)[0]
        if log_file_folder:
            os.makedirs(log_file_folder, exist_ok=True)
        fh = logging.FileHandler(log_file, 'w', encoding='utf-8')
        fh.setFormatter(logging.Formatter(fmt, datefmt))
        logger.addHandler(fh)
    logger.setLevel(log_level)
    return logger

fmt = '%(asctime)s %(filename)s %(levelname)s L %(lineno)d: %(message)s'
datefmt = '%Y-%m-%d %H:%M'
logger = get_logger()
# TODO esm model has limit on max seq len with batch size 5, at least 3000 is ok, 2700 is wrong.
MAX_PROTEIN_LEN = 2700


class ProtTransEmbedding(object):
    """ TODO just default t5_xl_half, use_fp16_in_pretrained_embedding is not fully realized.

    embedding length is 1024
    t5_xl_half is 3~5% better than bert_bfd according to info in the web below
    https://github.com/Rostlab/ProtTrans
    """
    def __init__(self, gpu_id=0, use_fp16_in_pretrained_embedding=True, mode_type='t5_xl_half') -> None:
        logger.info('Use ProtTransEmbedding')
        self.device = get_device(gpu_id)
        if mode_type == 't5_xl_half':
            self.model_root = '/mnt/sda/models/Rostlab/prot_t5_xl_half_uniref50-enc'
            self.tokenizer = T5Tokenizer.from_pretrained(self.model_root, do_lower_case=False)
            self.model = T5EncoderModel.from_pretrained(self.model_root).to(self.device)
            # only GPUs support half-precision currently; if you want to run on CPU use full-precision (not recommended, much slower)
            self.model.full() if self.device=='cpu' else self.model.half()
            self.model = self.model.eval()
            self.start_i = 0
            self.end_extra_i = 0
        elif mode_type == 'bert':
            self.model_root = '/mnt/sda/models/Rostlab/prot_bert_bfd'
            self.tokenizer = BertTokenizer.from_pretrained(self.model_root, do_lower_case=False)
            self.model = BertModel.from_pretrained(self.model_root).to(self.device)
            self.model = self.model.eval()
            self.start_i = 1
            self.end_extra_i = 1
        self.mode_type = mode_type

    def cal_pair_embeddings_and_convert2cpu(self, pep_seqs, prot_seqs, categories):
        """
        Returns:
            [prot_token_representations, pep_token_representations, category]
            Note: prot is the first
        Only consider t5_xl_half
        """
        dataset = []
        pep_token_representations = self.cal_embedding(pep_seqs)
        prot_token_representations = self.cal_embedding(prot_seqs)
        for i, (pep_seq, prot_seq, category) in enumerate(zip(pep_seqs, prot_seqs, categories)):
            dataset.append([
                # if consider both t5_xl_half and bert
                # pep_token_representations[i, self.start_i:len(pep_seq)+self.end_extra_i].cpu(),
                prot_token_representations[i, :len(prot_seq)].cpu(),
                pep_token_representations[i, :len(pep_seq)].cpu(),
                category
            ])
        return dataset

    def cal_embedding(self, seqs):
        """ Returns: a tensor on gpu"""
        _seqs = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in seqs]
        if self.mode_type == 'bert':
            ids = self.tokenizer(_seqs, return_tensors='pt')
        elif self.mode_type == 't5_xl_half':
            ids = self.tokenizer.batch_encode_plus(_seqs, add_special_tokens=True, padding="longest")
        input_ids = torch.LongTensor(ids['input_ids']).to(self.device)
        attention_mask = torch.LongTensor(ids['attention_mask']).to(self.device)
        with torch.no_grad():
            embedding_repr = self.model(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state
        return embedding_repr

    def cal_and_pad_embedding(self, seqs):
        """ Only consider t5_xl_half """
        token_representations = self.cal_embedding(seqs)
        tensors = []
        for i, seq in enumerate(seqs):
            # if consider both t5_xl_half and bert
            # tensors.append(token_representations[i, self.start_i:len(pep_seq)+self.end_extra_i])
            tensors.append(token_representations[i, :len(seq)])
        padded_tensor = pad_sequence(tensors, batch_first=True)
        return padded_tensor


class EsmEmbedding(object):
    """
    For 24GB GPU
    seq length 5000, batch size 1 is ok for esm model, but cannot be 5
    seq length 3000, batch size can be 5.
    """
    def __init__(self, gpu_id=0, use_fp16_in_pretrained_embedding=False, max_len=None) -> None:
        logger.info('Use EsmEmbedding')
        self.device = get_device(gpu_id)
        self.batch_tokenizer, self.model, self.esm_layer_num = get_esm_model(max_len=max_len)
        self.model.to(self.device)
        self.enable_fp16 = use_fp16_in_pretrained_embedding

    def cal_pair_embeddings_and_convert2cpu(self, pep_seqs, prot_seqs, categories):
        """
        Returns:
            [prot_token_representations, pep_token_representations, category]
            Note: prot is the first
        """
        dataset = []
        pep_token_representations = self.cal_embedding(pep_seqs)
        prot_token_representations = self.cal_embedding(prot_seqs)
        for i, (pep_seq, prot_seq, category) in enumerate(zip(pep_seqs, prot_seqs, categories)):
            start_i = 1
            pep_end_i = len(pep_seq) + 1
            prot_end_i = len(prot_seq) + 1
            dataset.append([
                prot_token_representations[i, start_i: prot_end_i].cpu(),
                pep_token_representations[i, start_i: pep_end_i].cpu(),
                category,
            ])
        return dataset

    def cal_embedding(self, seqs):
        """ keeps gpu if gpu is available """
        batch_tokens = self.batch_tokenizer(seqs)
        batch_tokens = batch_tokens.to(self.device)
        with torch.no_grad():
            results = self.model(batch_tokens, repr_layers=[self.esm_layer_num], return_contacts=False)
        token_representations = results["representations"][self.esm_layer_num]
        if self.enable_fp16:
            token_representations = token_representations.to(torch.float16)
        return token_representations

    def cal_and_pad_embedding(self, seqs):
        """ After cal_embedding, the for each seq, each pad_idx aa has different 1d values, and so remove padded idx
        embeddings at end of seq """
        token_representations = self.cal_embedding(seqs)
        tensors = []
        for i, seq in enumerate(seqs):
            start_i = 1
            prot_end_i = len(seq) + 1
            tensors.append(token_representations[i, start_i: prot_end_i])
        padded_tensor = pad_sequence(tensors, batch_first=True)
        return padded_tensor


class ESMBatchTokenizer(object):
    """Callable to convert an unprocessed (labels + strings) batch to a
    processed (labels + tensor) batch.
    """

    def __init__(self, alphabet, truncation_seq_length: int = None, max_len=50):
        self.alphabet = alphabet
        self.truncation_seq_length = truncation_seq_length
        self.max_len = max_len

    def __call__(self, seq_str_list: list[str]):
        batch_size = len(seq_str_list)
        seq_encoded_list = [self.alphabet.encode(seq_str) for seq_str in seq_str_list]
        if self.truncation_seq_length:
            seq_encoded_list = [seq_str[:self.truncation_seq_length] for seq_str in seq_encoded_list]
        if isinstance(self.max_len, int) and self.max_len > 0:
            max_len = self.max_len
        else:
            max_len = max(len(seq_encoded) for seq_encoded in seq_encoded_list)
        tokens = torch.empty(
            (
                batch_size,
                max_len + int(self.alphabet.prepend_bos) + int(self.alphabet.append_eos),
            ),
            dtype=torch.int64,
        )
        tokens.fill_(self.alphabet.padding_idx)

        for i, seq_encoded in enumerate(seq_encoded_list):
            if self.alphabet.prepend_bos:
                tokens[i, 0] = self.alphabet.cls_idx
            seq = torch.tensor(seq_encoded, dtype=torch.int64)
            tokens[
                i,
                int(self.alphabet.prepend_bos) : len(seq_encoded)
                + int(self.alphabet.prepend_bos),
            ] = seq
            if self.alphabet.append_eos:
                tokens[i, len(seq_encoded) + int(self.alphabet.prepend_bos)] = self.alphabet.eos_idx

        return tokens


@cache
def get_esm_model(model_level='650m', max_len=None):
    if model_level == '650m':
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        esm_layer_num = 33
    elif model_level == '3b':
        model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
        esm_layer_num = 36
    batch_tokenizer = ESMBatchTokenizer(alphabet, max_len=max_len)
    model.eval()
    return batch_tokenizer, model, esm_layer_num


if __name__ == "__main__":
    from icecream import ic
    ic.configureOutput(includeContext=True, argToStringFunction=lambda _: str(_))
    import os, sys
    sys.path.append(os.path.abspath('.'))
    from utils.file_util import FileUtil


    seqs = [
        'ADE',
        'ADEAA',
    ]
    use_esm = 1
    if use_esm:
        embedder = EsmEmbedding(gpu_id=1)
        ic(embedder.batch_tokenizer.alphabet.prepend_bos)
        ic(embedder.batch_tokenizer.alphabet.append_eos)
        ic(embedder.batch_tokenizer.alphabet.padding_idx)
        ic(embedder.batch_tokenizer.alphabet.cls_idx)
        ic(embedder.batch_tokenizer.alphabet.eos_idx)
    else:
        embedder = ProtTransEmbedding(gpu_id=1, mode_type='bert')  # bert t5_xl_half

    for length in range(2700, 2800, 100):
        # length = 3000
        logger.info(f'length {length}')
        long_seqs = [''.join(['A'] * length)] * 5
        # ic(len(long_seqs), long_seqs)
        # long_seqs = FileUtil.read_raw_text('1.txt')[:16]
        # lens = [len(seq) for seq in long_seqs]
        # ic(max(lens))
        result = embedder.cal_embedding(long_seqs)
        logger.info(f'type(result) {type(result)}')
        logger.info(f' {result.shape}')
        logger.info(f' {result[0]}')