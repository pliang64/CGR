import csv
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from PIL import Image
from Bio import SeqIO
from sklearn.preprocessing import RobustScaler, MinMaxScaler
from matplotlib.colors import LogNorm

def parse_sequences(filename):
    sequences = []
    with open(filename) as inputfile:
        for record in SeqIO.parse(inputfile, "fasta"):
            sequence_name = record.id
            sequence = str(record.seq)
            sequences.append((sequence_name, sequence))
    return sequences

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

def generate_images(filepath, size=512, output_folder='./images/',
    logs_folder='./logs/', int_matrix=False, affix='', show_plot=False, robust_scaler=False, no_norm=False, cmap='gray'):

    # get file name and folder path
    file, folder = os.path.basename(filepath), os.path.dirname(filepath)
    output_files = []

    # # Load file and plot coordinates
    # dna = parse_sequence(folder+"/"+file)
    # coords = dna_cgr(dna)
    
    sequences = parse_sequences(folder+"/"+file)

    if robust_scaler:
        rbst = '_robust'
    else:
        rbst = ''
        
    for sequence_name, sequence in sequences:
        output_file = output_folder + sequence_name + affix + '_' + str(size) + rbst + "_" + cmap
        output_files.append(output_file) # Append output even if it exists
        # Check if output_file exists, skips if true
        if os.path.exists(output_file + '_cgr.png'):
            #print(f"{output_file} already exists. Skipping.")
            continue
        
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
            # plt.imshow(img_inv, cmap=cmap, norm=norm)
            plt.imshow(img_inv, cmap=cmap)
            # plt.imshow(img_norm, cmap=cmap)
            plt.show()
            # plt.savefig(output_file+'_colormap.png')
            plt.close()
            
        
        Image.fromarray(img_inv).save(output_file + '_cgr.png')
        # Image.fromarray(img_norm).save(output_file + '_cgr.png')
        
        # end timer
        end_time = time.time()
        elapsed_time = end_time - start_time

        print(f"{output_file}: Elapsed time: {elapsed_time:.2f} seconds")
        # with open(logs_folder+'processing_time.txt', 'a') as f:
        #     f.write(f"{file}: Elapsed time: {elapsed_time:.2f} seconds\n")
    #return img_rot
    return output_files