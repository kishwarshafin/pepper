import shutil
import pickle
from os import listdir, remove
from os.path import isfile, join
"""
EXAMPLE USAGE:
bed_directory_path = "/Users/saureous/data/bed_alleles_copy"
output_file_path = "output/bed_output/concatenated.bed"

file_manager = FileManager()
file_paths = file_manager.get_file_paths_from_directory(directory_path=bed_directory_path)
file_manager.concatenate_files(file_paths=file_paths, output_file_path=output_file_path)
file_manager.delete_files(file_paths=file_paths)
"""


class FileManager:
    """
    Does simple file operations like concatenation, fetching a list of paths for files in a directory, deleting files
    """
    @staticmethod
    def concatenate_files(file_paths, output_file_path):
        """
        Concatenate files given in list of file paths to a single file
        :param file_paths: List of file path
        :param output_file_path: Output file path name
        :return:
        """
        with open(output_file_path, 'wb') as out_file:
            for file_path in file_paths:
                with open(file_path, 'rb') as in_file:
                    # 100MB per writing chunk to avoid reading big file into memory.
                    shutil.copyfileobj(in_file, out_file, 1024*1024*100)

    @staticmethod
    def merge_dictionaries(file_paths, output_file_path):
        """
        Concatenate dictionaries given in list of file paths to a single dictionary
        :param file_paths: List of file path
        :param output_file_path: Output file path name
        :return:
        """
        merged_dict = {}

        for file_path in file_paths:
            dictionary = pickle.load(open(file_path, "rb"))
            z = {**merged_dict, **dictionary}
            merged_dict = z

        # save the merged file
        pickle.dump(merged_dict, open(output_file_path, "wb"))

    @staticmethod
    def get_file_paths_from_directory(directory_path):
        """
        Returns all paths of files given a directory path
        :param directory_path: Path to the directory
        :return: A list of paths of files
        """
        file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))]
        return file_paths

    @staticmethod
    def delete_files(file_paths):
        """
        Deletes files given in file paths
        :param file_paths: List of file paths
        :return:
        """
        for file_path in file_paths:
            remove(file_path)

