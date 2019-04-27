from torch.utils.data import Dataset
import torchvision.transforms as transforms
import h5py


class SequenceDataset(Dataset):
    """
    Arguments:
        A HDF5 file path
    """
    def __init__(self, hdf5_file_path):
        self.transform = transforms.Compose([transforms.ToTensor()])
        file_image_pair = []

        with h5py.File(hdf5_file_path, 'r') as hdf5_file:
            image_names = list(hdf5_file['summaries'].keys())

        for image_name in image_names:
            file_image_pair.append((hdf5_file_path, image_name))

        self.all_images = file_image_pair

    def __getitem__(self, index):
        # load the image
        hdf5_filepath, image_name = self.all_images[index]

        with h5py.File(hdf5_filepath, 'r') as hdf5_file:
            images = hdf5_file['summaries'][image_name]['image'][()]
            position = hdf5_file['summaries'][image_name]['position'][()]
            index = hdf5_file['summaries'][image_name]['index'][()]
            contig = hdf5_file['summaries'][image_name]['contig'][()]
            chunk_id = hdf5_file['summaries'][image_name]['chunk_id'][()]
            contig_start = hdf5_file['summaries'][image_name]['region_start'][()]
            contig_end = hdf5_file['summaries'][image_name]['region_end'][()]
            hp_tag = hdf5_file['summaries'][image_name]['hp'][()]

        return contig, contig_start, contig_end, chunk_id, images, position, index, hp_tag

    def __len__(self):
        return len(self.all_images)
