import h5py
import yaml
import numpy as np
from pepper_variant.modules.python.Options import ImageSizeOptions


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

    def write_summary(self, summary_name, images, positions, base_labels, type_labels, contig_name, region_start, region_end):
        if 'summaries' not in self.meta:
            self.meta['summaries'] = set()

        if summary_name not in self.meta['summaries']:
            self.meta['summaries'].add(summary_name)
            self.file_handler.create_group('{}/{}'.format(self._summary_path_, summary_name))
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'region_start')] = region_start
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'region_end')] = region_end
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, 'contig')] = contig_name

        # create the group
        image_size = (len(images), len(images[0]), len(images[0][0]))
        position_size = (len(positions), 1)
        base_size = (len(base_labels), 1)
        type_size = (len(type_labels), 1)

        self.file_handler['{}/{}'.format(self._summary_path_, summary_name)].create_dataset("images", image_size, dtype=np.int, compression="gzip", data=images)
        self.file_handler['{}/{}'.format(self._summary_path_, summary_name)].create_dataset("positions", position_size, dtype=np.int32, data=positions)
        self.file_handler['{}/{}'.format(self._summary_path_, summary_name)].create_dataset("base_labels", base_size, dtype=np.uint8, data=base_labels)
        self.file_handler['{}/{}'.format(self._summary_path_, summary_name)].create_dataset("type_labels", type_size, dtype=np.uint8, data=type_labels)

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
