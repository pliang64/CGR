import csv
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from PIL import Image
from Bio import SeqIO
from sklearn.preprocessing import RobustScaler, MinMaxScaler
from matplotlib.colors import LogNorm
import argparse

def get_concatenated_sequence(species_name, sample_name, folder_path):
    concatenated_sequence = ''

    # Read files and concatenate sequences using BioPython
    for chr_num in range(1, 20):
        file_name = f"{species_name}.{chr_num}.fa"
        file_path = os.path.join(folder_path, file_name)

        if not os.path.exists(file_path):
            print(f"File {file_name} does not exist in the folder {folder_path}.")
            return None

        for record in SeqIO.parse(file_path, "fasta"):
            if record.id == sample_name:
                concatenated_sequence += str(record.seq)
                break

    return concatenated_sequence

def dna_cgr(seq):
    N = len(seq)
    dataPoints = np.zeros((2, N+1))
    dataPoints[:, 0] = np.array([0.5, 0.5])
    for i in range(1, N+1):
        if seq[i-1] == 'A':
            corner = np.array([0, 0])
        elif seq[i-1] == 'C':
            corner = np.array([0, 1])
        elif seq[i-1] == 'G':
            corner = np.array([1, 1])
        elif seq[i-1] == 'T':
            corner = np.array([1, 0])
        else:
            corner = np.array([0.5, 0.5])
        dataPoints[:, i] = 0.5*(dataPoints[:, i-1] + corner)
    return(dataPoints)

def generate_images(sequence_name, sequence, size=512, output_folder='./images/',
    logs_folder='./logs/', int_matrix=False, affix='', show_plot=False, robust_scaler=False, no_norm=False, cmap='gray'):


    if robust_scaler:
        rbst = '_robust'
    else:
        rbst = ''
 
    output_file = output_folder + sequence_name + affix + '_' + str(size) + rbst
    if robust_scaler:
        output_file += "_robust"
    elif not no_norm:
        output_file += "_lognorm"
        
    # Check if output_file exists, skips if true
    if os.path.exists(output_file + '_cgr.png'):
        print(f"{output_file} already exists. Skipping.")
        return output_file

    coords = dna_cgr(sequence)
    print(f"Sequence: {sequence_name}")


    # start timer
    start_time = time.time()

    # Create a 2D grayscale image with size: sizexsize
    #img = np.zeros((size, size), dtype=np.uint16)  # Using uint16 data type, will convert after normalizing
    img = np.zeros((size, size), dtype=np.uint32)  # Using uint8 data type, will convert after normalizing

    max_value = 0;
    # Set pixels in the image based on the coordinates
    for i in range(len(coords[0])):
        x = int(coords[0][i] * size - 1)
        y = int(coords[1][i] * size - 1)
        img[x, y] += 1
        if img[x, y] > max_value:
            max_value = img[x, y]
        #print("x:",x,", y:",y)
    print(f"Max Value: '{max_value}'")

    # Rotate the image counterclockwise by 90 degrees so that [0,0] corresponds with bottom left
    img_rot = np.rot90(img)

    # Apply the RobustScaler

    if robust_scaler:
        robust_scaler = RobustScaler()
        robust_data = robust_scaler.fit_transform(img_rot)
        # Apply the MinMaxScaler to the scaled values from RobustScaler
        min_max_scaler = MinMaxScaler(feature_range=(0, 1))
        img_norm = (min_max_scaler.fit_transform(robust_data)*255).astype(np.uint8)
    elif no_norm:
        img_norm = img_rot.astype(np.uint8)
    else:
        # img_norm = (img_rot / np.max(img_rot) * 255).astype(np.uint8)
        img_norm = np.log1p(img_rot) 
        img_norm = (img_norm - np.min(img_norm)) / (np.max(img_norm) - np.min(img_norm)) * 255
        img_norm = img_norm.astype(np.uint8)

        # img_norm = img_rot
        # Create a normalized colormap object to map values to colors

    if int_matrix:
        np.savetxt(output_file + "_int_matrix.txt", img_rot, delimiter=",", fmt="%d")

    # Invert the image
    img_inv = 255 - img_norm


    # Save and display the CGR representation
    if show_plot:
        plt.imshow(img_inv, cmap=cmap)
        plt.show()
        plt.close()


    Image.fromarray(img_inv).save(output_file + '_cgr.png')

    # end timer
    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"{output_file}: Elapsed time: {elapsed_time:.2f} seconds")
    return output_file

def process_all_samples(species_name, folder_path, output_path):
    sample_names = []

    # Extract sample names from the first chromosome file
    first_chr_file = f"{species_name}.1.fa"
    first_chr_file_path = os.path.join(folder_path, first_chr_file)

    if not os.path.exists(first_chr_file_path):
        print(f"File {first_chr_file} does not exist in the folder {folder_path}.")
        return

    for record in SeqIO.parse(first_chr_file_path, "fasta"):
        sample_names.append(record.id)

    # print(sample_names)
    
    # Perform get_concatenated_sequence for each sample_name
    for sample_name in sample_names:
        concat_sequence = get_concatenated_sequence(species_name, sample_name, folder_path)
        if concat_sequence:
            print(f"Concatenated sequence for {sample_name}.")
            generate_images(sample_name, concat_sequence, size=224, output_folder=output_path, show_plot=True, robust_scaler=True)
        else:
            print(f"Could not get the concatenated sequence for {sample_name}.")

def main(species_name, folder_path, output_path):
    output_path = os.path.join(output_path,species_name)
    os.makedirs(output_path, exist_ok=True)
    process_all_samples(species_name, folder_path, output_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some arguments.')
    parser.add_argument('--species_name', type=str, required=True, help='Name of the species.')
    parser.add_argument('--folder_path', type=str, required=True, help='Path to the folder containing FASTA files.')
    parser.add_argument('--output_path', type=str, required=True, help='Path to save the output CGR files.')

    args = parser.parse_args()

    main(args.species_name, args.folder_path, args.output_path)

