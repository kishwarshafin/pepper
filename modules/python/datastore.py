import h5py
import yaml
import numpy as np
from modules.python.Options import ImageSizeOptions


class DataStore(object):
    """Class to read/write to a HELEN's file"""
    _summary_path_ = 'summaries'
    _groups_ = ('image', 'position', 'index', 'label')

    def __init__(self, filename, mode='r'):
        self.filename = filename
        self.mode = mode

        self._sample_keys = set()
        self.file_handler = h5py.File(self.filename, self.mode)

        self._meta = None

    def __enter__(self):
        self.file_handler = None

        return self

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

    def write_train_summary(self, contig_name, image, label, position, index, summary_name):
        self.file_handler['{}/{}/{}/{}'.format(self._summary_path_, contig_name, summary_name, 'image')] = np.array(image, dtype=np.float)
        self.file_handler['{}/{}/{}/{}'.format(self._summary_path_, contig_name, summary_name, 'label')] = np.array(label, dtype=np.uint8)
        self.file_handler['{}/{}/{}/{}'.format(self._summary_path_, contig_name, summary_name, 'position')] = position
        self.file_handler['{}/{}/{}/{}'.format(self._summary_path_, contig_name, summary_name, 'index')] = np.array(index, dtype=np.uint8)
