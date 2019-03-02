import math
import argparse
import torch
import numpy as np
import h5py
from torchvision import transforms


def get_strand_color(is_rev):
    """
    Get color for forward and reverse reads
    :param is_rev: True if read is reversed
    :return:
    """
    is_rev = int(math.ceil(is_rev))
    if is_rev == 254:
        return 'R'
    if is_rev == 240:
        return '1'
    elif is_rev == 70:
        return '0'
    else:
        return ' '


def get_alt_type(alt_type_color):
    """
    Get color for forward and reverse reads
    :param is_rev: True if read is reversed
    :return:
    """
    alt_type_color = int(math.ceil(alt_type_color))
    if alt_type_color == 0:
        return ' '
    elif alt_type_color == 5:
        return '0'
    elif alt_type_color == 240:
        return '1'
    elif alt_type_color == 125:
        return '2'
    elif alt_type_color == 254:
        return 'R'


def get_base_from_color(base_color):
    color = int(math.ceil(base_color))
    global_base_color_reverse = {200: 'A', 50: 'C', 150: 'G', 100: 'T', 10: '.', 250: '*'}
    # {{'C', 50}, {'T', 100}, {'G', 150}, {'A', 200}, {'*', 250}, {'.', 10}, {'N', 10}};
    if color in global_base_color_reverse:
        return global_base_color_reverse[color]
    else:
        return ' '
    # 'A': 25.0, 'C': 75.0, 'G': 125.0, 'T': 175.0, '*': 225.0


def get_quality_by_color(quality):
    """
    Get a color spectrum given mapping quality
    :param map_quality: value of mapping quality
    :return:
    """
    quality = int(math.ceil(quality))
    color = math.floor(((quality / 254) * 9))
    if color == 0:
        return ' '
    return str(color)


def get_mismatch_or_alt_color(alt_type_color):
    """
    Get color for forward and reverse reads
    :param is_rev: True if read is reversed
    :return:
    """
    alt_type_color = int(math.ceil(alt_type_color))
    if alt_type_color == 0:
        return ' '
    elif alt_type_color == 50:
        return '0'
    elif alt_type_color == 250:
        return '1'


def analyze_tensor(image):
    # base_color, base_quality_color, map_qual_color, strand_color, alt_color
    img_c, img_w, img_h = image.size()
    image = np.array(image.data * 254)
    img_h = 100
    print("BASE CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_base_from_color(image[0][j][i]), end='')
        print()

    print("BASE QUALITY CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_quality_by_color(image[1][j][i]), end='')
        print()
    print("MAPPING QUALITY CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_quality_by_color(image[2][j][i]), end='')
        print()
    print("STRAND DIRECTION CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_strand_color(image[3][j][i]), end='')
        print()
    print("MISMATCH CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_mismatch_or_alt_color(image[4][j][i]), end='')
        print()
    print("ALT1 CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_mismatch_or_alt_color(image[5][j][i]), end='')
        print()
    print("ALT2 CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_mismatch_or_alt_color(image[6][j][i]), end='')
        print()


def save_base_quality_array(image):
    img_c, img_w, img_h = image.size()
    image = np.array(image.data * 254)
    img_h = 100
    entire_image = []
    for i in range(img_h):
        image_row = []
        for j in range(img_w):
            image_row.append([0, 0, image[1][j][i], 255])
        entire_image.append(image_row)
        print()
    entire_image = np.array(entire_image)
    from scipy import misc
    misc.imsave("base_quality_tensor" + ".png", entire_image, format="PNG")


def save_map_quality_array(image):
    img_c, img_w, img_h = image.size()
    image = np.array(image.data * 254)
    img_h = 100
    entire_image = []
    for i in range(img_h):
        image_row = []
        for j in range(img_w):
            image_row.append([0, image[2][j][i], 0, 255])
        entire_image.append(image_row)
        print()
    entire_image = np.array(entire_image)
    from scipy import misc
    misc.imsave("map_quality_tensor" + ".png", entire_image, format="PNG")


def save_strand_array(image):
    img_c, img_w, img_h = image.size()
    image = np.array(image.data * 254)
    img_h = 100
    entire_image = []
    for i in range(img_h):
        image_row = []
        for j in range(img_w):
            image_row.append([image[3][j][i], 0, 0, 255])
        entire_image.append(image_row)
        print()
    entire_image = np.array(entire_image)
    from scipy import misc
    misc.imsave("strand_color_tensor" + ".png", entire_image, format="PNG")


def save_mismatch(image):
    img_c, img_w, img_h = image.size()
    image = np.array(image.data * 254)
    img_h = 100
    entire_image = []
    for i in range(img_h):
        image_row = []
        for j in range(img_w):
            image_row.append([0, image[4][j][i], 0, 255])
        entire_image.append(image_row)
        print()
    entire_image = np.array(entire_image)
    from scipy import misc
    misc.imsave("mismatch_tensor" + ".png", entire_image, format="PNG")


def save_alt1(image):
    img_c, img_w, img_h = image.size()
    image = np.array(image.data * 254)
    img_h = 100
    entire_image = []
    for i in range(img_h):
        image_row = []
        for j in range(img_w):
            image_row.append([image[5][j][i], 0, 0, 255])
        entire_image.append(image_row)
        print()
    entire_image = np.array(entire_image)
    from scipy import misc
    misc.imsave("alt1_freq_tensor" + ".png", entire_image, format="PNG")


def save_alt2(image):
    img_c, img_w, img_h = image.size()
    image = np.array(image.data * 254)
    img_h = 100
    entire_image = []
    for i in range(img_h):
        image_row = []
        for j in range(img_w):
            image_row.append([0, 0, image[6][j][i], 255])
        entire_image.append(image_row)
        print()
    entire_image = np.array(entire_image)
    from scipy import misc
    misc.imsave("alt2_freq_tensor" + ".png", entire_image, format="PNG")


def save_base_array(image):
    img_c, img_w, img_h = image.size()
    image = np.array(image.data * 254)
    img_h = 100
    entire_image = []
    for i in range(img_h):
        image_row = []
        for j in range(img_w):
            if image[0][j][i] != 0:
                print(get_base_from_color(image[0][j][i]), end='')
                if get_base_from_color(image[0][j][i]) == ' ':
                    image_row.append([255, 255, 255, 255])
                # elif get_base_from_color(image[0][j][i]) == get_base_from_color(image[0][j][0]) and i > 0:
                #     image_row.append([255, 255, 255, 255])
                elif get_base_from_color(image[0][j][i]) == 'A':
                    image_row.append([0, 0, 255, 255])
                elif get_base_from_color(image[0][j][i]) == 'C':
                    image_row.append([255, 0, 0, 255])
                elif get_base_from_color(image[0][j][i]) == 'G':
                    image_row.append([0, 255, 0, 255])
                elif get_base_from_color(image[0][j][i]) == 'T':
                    image_row.append([255, 255, 0, 255])
                else:
                    # purple
                    image_row.append([160, 32, 240, 255])
            else:
                print(' ', end='')
                image_row.append([0, 0, 0, 255])
        entire_image.append(image_row)
        print()
    entire_image = np.array(entire_image)
    from scipy import misc
    misc.imsave("base_tensor" + ".png", entire_image, format="PNG")


def tensor_to_image(image):
    # base_color, base_quality_color, map_qual_color, strand_color, alt_color
    print("BASE CHANNEL:")
    save_base_array(image)
    save_base_quality_array(image)
    save_map_quality_array(image)
    save_strand_array(image)
    save_mismatch(image)
    save_alt1(image)
    save_alt2(image)


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tensor_file",
        type=str,
        required=True,
        help="H5PY file path"
    )
    parser.add_argument(
        "--index",
        type=int,
        required=True,
        help="Index of image."
    )
    FLAGS, unparsed = parser.parse_known_args()

    if FLAGS.tensor_file:
        hdf5_image = FLAGS.tensor_file
        hdf5_index = FLAGS.index
        hdf5_file = h5py.File(hdf5_image, 'r')

        image_dataset = hdf5_file['images']
        image = np.array(image_dataset[hdf5_index], dtype=np.uint8)
        transform = transforms.Compose([transforms.ToTensor()])
        image = transform(image)
        image = image.transpose(1, 2)

        label_dataset = hdf5_file['labels']
        label = np.array(label_dataset[hdf5_index], dtype=np.long)
        for l in label:
            print(l, end='')
        print()
        # analyze_tensor(image)

        tensor_to_image(image)

