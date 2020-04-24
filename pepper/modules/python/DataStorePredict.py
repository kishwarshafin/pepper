import h5py
import yaml
import numpy as np


class DataStore(object):
    """Class to read/write to a FRIDAY's file"""
    _prediction_path_ = 'predictions'
    _groups_ = ('position', 'index', 'bases', 'rles')

    def __init__(self, filename, mode='r'):
        self.filename = filename
        self.mode = mode

        self._sample_keys = set()
        self.file_handler = h5py.File(self.filename, self.mode)

        self._meta = None

    def __exit__(self, *args):
        if self.mode != 'r' and self._meta is not None:
            self._write_metadata(self.meta)
        self.file_handler.close()

    def _write_metadata(self, data):
        """Save a data structure to file within a yml str."""
        for group, d in data.items():
            if group in self.file_handler:
                del self.file_handler[group]
            self.file_handler[group] = yaml.dump(d)

    def _load_metadata(self, groups=None):
        """Load meta data"""
        if groups is None:
            groups = self._groups_
        return {g: yaml.load(self.file_handler[g][()]) for g in groups if g in self.file_handler}

    @property
    def meta(self):
        if self._meta is None:
            self._meta = self._load_metadata()
        return self._meta

    def update_meta(self, meta):
        """Update metadata"""
        self._meta = self.meta
        self._meta.update(meta)

    def write_prediction(self, contig, contig_start, contig_end, chunk_id, position, index, predicted_bases, phred_score):
        chunk_name_prefix = str(contig) + "-" + str(contig_start.item()) + "-" + str(contig_end.item())
        chunk_name_suffix = str(chunk_id.item())

        name = contig + chunk_name_prefix + chunk_name_suffix

        if 'predictions' not in self.meta:
            self.meta['predictions'] = set()
        if 'predictions_contig' not in self.meta:
            self.meta['predictions_contig'] = set()

        if chunk_name_prefix not in self.meta['predictions_contig']:
            self.meta['predictions_contig'].add(chunk_name_prefix)
            self.file_handler['{}/{}/{}/{}'.format(self._prediction_path_, contig, chunk_name_prefix, 'contig_start')] \
                = contig_start.item()
            self.file_handler['{}/{}/{}/{}'.format(self._prediction_path_, contig, chunk_name_prefix, 'contig_end')] \
                = contig_end.item()

        if name not in self.meta['predictions']:
            self.meta['predictions'].add(name)
            self.file_handler['{}/{}/{}/{}/{}'.format(self._prediction_path_, contig, chunk_name_prefix,
                                                      chunk_name_suffix, 'position')] = position
            self.file_handler['{}/{}/{}/{}/{}'.format(self._prediction_path_, contig, chunk_name_prefix,
                                                      chunk_name_suffix, 'index')] = index
            self.file_handler['{}/{}/{}/{}/{}'.format(self._prediction_path_, contig, chunk_name_prefix,
                                                      chunk_name_suffix, 'bases')] = predicted_bases.astype(np.uint8)
            self.file_handler['{}/{}/{}/{}/{}'.format(self._prediction_path_, contig, chunk_name_prefix,
                                                      chunk_name_suffix, 'phred_score')] = phred_score.astype(np.uint8)


