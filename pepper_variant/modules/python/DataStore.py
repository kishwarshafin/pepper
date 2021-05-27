import h5py
import yaml
import numpy as np


class DataStore(object):
    """Class to read/write to a HELEN's file"""
    _summary_path_ = 'summaries'
    _groups_ = ('image', 'position', 'index', 'label')

    def __init__(self, filename, mode='r'):
        self.filename = filename
        self.mode = mode

        self._sample_keys = set()
        self.file_handler = None

        self._meta = None

    def __enter__(self):
        self.file_handler = h5py.File(self.filename, self.mode)
        return self

    def __exit__(self, *args):
        # if self.mode != 'r' and self._meta is not None:
        #     self._write_metadata(self.meta)
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

    def write_summary(self, summary_name, candidate):
        if 'summaries' not in self.meta:
            self.meta['summaries'] = set()

        if summary_name not in self.meta['summaries']:
            self.meta['summaries'].add(summary_name)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'image')] = np.array(candidate.image_matrix, dtype=np.uint8)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'base_label')] = np.array(candidate.base_label, dtype=np.uint8)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'type_label')] = np.array(candidate.type_label, dtype=np.uint8)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'position')] = np.array(candidate.position, dtype=np.int32)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'region_start')] = np.array(candidate.region_start, dtype=np.int32)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'region_stop')] = np.array(candidate.region_stop, dtype=np.int32)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'contig')] = candidate.contig

    def write_summary_hp(self, region, image_hp1, image_hp2, label_hp1, label_hp2, position, index, chunk_id, summary_name):
        contig_name, region_start, region_end = region
        if 'summaries' not in self.meta:
            self.meta['summaries'] = set()

        if summary_name not in self.meta['summaries']:
            self.meta['summaries'].add(summary_name)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'image_hp1')] = np.array(image_hp1, dtype=np.uint8)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'image_hp2')] = np.array(image_hp2, dtype=np.uint8)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'label_hp1')] = np.array(label_hp1, dtype=np.uint8)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'label_hp2')] = np.array(label_hp2, dtype=np.uint8)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'position')] = np.array(position, dtype=np.int32)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'index')] = np.array(index, dtype=np.int32)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'contig')] = contig_name
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'region_start')] = region_start
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'region_end')] = region_end
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'chunk_id')] = chunk_id
