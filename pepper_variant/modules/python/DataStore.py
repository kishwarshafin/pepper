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

    def write_summary(self, summary_name, contigs, positions, depths, all_candidates, all_candidate_frequency, all_images, all_base_labels, all_type_label, train_mode):
        if 'summaries' not in self.meta:
            self.meta['summaries'] = set()

        # create the group

        dt_candidates = h5py.special_dtype(vlen=str)
        if summary_name not in self.meta['summaries']:
            self.meta['summaries'].add(summary_name)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, "contigs")] = np.array(contigs, dtype='S')
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, "positions")] = np.array(positions, dtype=np.int32)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, "depths")] = np.array(depths, dtype=np.uint8)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, "candidates")] = np.array(all_candidates, dtype=dt_candidates)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, "candidate_frequency")] = np.array(all_candidate_frequency, dtype=np.uint8)
            self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, "images")] = np.array(all_images, dtype=np.int8)
            if train_mode:
                self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, "base_labels")] = np.array(all_base_labels, dtype=np.uint8)
                self.file_handler['{}/{}/{}'.format(self._summary_path_, summary_name, "type_label")] = np.array(all_type_label, dtype=np.uint8)
        # self.file_handler['{}'.format(self._summary_path_)].create_dataset("contigs", contig_size, dtype=str_dt, data=contigs)
        # self.file_handler['{}'.format(self._summary_path_)].create_dataset("positions", position_size, dtype=np.int32, data=positions)
        # self.file_handler['{}'.format(self._summary_path_)].create_dataset("depth", depth_size, dtype=np.int32, data=depths)
        # self.file_handler['{}'.format(self._summary_path_)].create_dataset("depth", candidate_frequency_size, dtype=np.int32, data=depths)
        #
        # self.file_handler['{}'.format(self._summary_path_)].create_dataset("images", image_size, dtype=np.int, compression="gzip", data=images)
        # self.file_handler['{}'.format(self._summary_path_)].create_dataset("images", image_size, dtype=np.int, compression="gzip", data=images)
        #
        # self.file_handler['{}'.format(self._summary_path_)].create_dataset("base_labels", base_size, dtype=np.uint8, data=base_labels)
        # self.file_handler['{}'.format(self._summary_path_)].create_dataset("type_labels", type_size, dtype=np.uint8, data=type_labels)

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
