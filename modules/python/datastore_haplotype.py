import h5py
import yaml
from collections import Counter, defaultdict, OrderedDict


class DataStore(object):
    """Class to read/write to a HELEN's file"""
    _prediction_path_ = 'predictions'
    _groups_ = ('predictions', 'position')

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

    def write_prediction(self, contig, contig_start, contig_end, chunk_id, hp_tag, position, index, predictions):
        contig_start = contig_start.item()
        contig_end = contig_end.item()
        chunk_name_prefix = str(contig) + "-" + str(contig_start) + "-" + str(contig_end)
        chunk_name_suffix = str(chunk_id.item())
        hp_tag = str(hp_tag.item())
        if 'chunk_name_prefix' not in self.meta:
            self.meta['chunk_name_prefix'] = set()

        if chunk_name_prefix not in self.meta['chunk_name_prefix']:
            self.meta['chunk_name_prefix'].add(chunk_name_prefix)
            self.file_handler['{}/{}/{}/{}'.format(self._prediction_path_, contig,
                                                   chunk_name_prefix, 'contig')] = contig
            self.file_handler['{}/{}/{}/{}'.format(self._prediction_path_, contig,
                                                   chunk_name_prefix, 'region_start')] = contig_start
            self.file_handler['{}/{}/{}/{}'.format(self._prediction_path_, contig,
                                                   chunk_name_prefix, 'region_end')] = contig_end

        self.file_handler['{}/{}/{}/{}/{}/{}'.format(self._prediction_path_, contig, chunk_name_prefix, hp_tag,
                                                     chunk_name_suffix, 'position')] = position
        self.file_handler['{}/{}/{}/{}/{}/{}'.format(self._prediction_path_, contig, chunk_name_prefix, hp_tag,
                                                     chunk_name_suffix, 'index')] = index
        self.file_handler['{}/{}/{}/{}/{}/{}'.format(self._prediction_path_, contig, chunk_name_prefix, hp_tag,
                                                     chunk_name_suffix, 'predictions')] = predictions
