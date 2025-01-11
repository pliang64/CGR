from cgr_by_chrom_v3 import generate_images
from process_folders_species import process_files
import os
import argparse

# Initialize argument parser
parser = argparse.ArgumentParser(description='Process specified input file for CGR generation.')
parser.add_argument('--input_file', type=str, help='Path to the input .fa file')
parser.add_argument('--out_folder', type=str, help='Path to the output folder')
parser.add_argument('--size', type=int, default=512, help='Size of the output image')
parser.add_argument('--norm', type=str, default='no_norm', choices=['robust_scaler', 'no_norm'], help='Normalization option')


# Parse the arguments
args = parser.parse_args()

# Assign the arguments to variables
input_file = args.input_file
out_folder = args.out_folder
size = args.size
norm = args.norm

# Other code logic stays the same
print(input_file)

# Check by chrom
filename = os.path.basename(input_file)
prefix = filename.split(".")[0]
#prefix = prefix.split("_")[0]
chr_num = filename.split(".")[1]
chr_num = "chr" + chr_num

# Set folder path
folder_path = os.path.join(out_folder, prefix, chr_num, "")
os.makedirs(folder_path, exist_ok=True)

# CGR
print(f"Processing '{filename}':")
print(folder_path)

# Determine the scaler
robust_scaler = True if norm == 'robust_scaler' else False

# # Generate images
# if robust_scaler:
#     generate_images(input_file, size=size, int_matrix=False, output_folder=folder_path, affix=("_" + chr_num), show_plot=False, robust_scaler=True)
# else:
#     generate_images(input_file, size=size, int_matrix=False, output_folder=folder_path, affix=("_" + chr_num), show_plot=False, no_norm=True)
    
# Generate images
if robust_scaler:
    generate_images(input_file, size=size, int_matrix=False, output_folder=folder_path, affix=("_" + chr_num), show_plot=False, robust_scaler=True, sequence_names_file="./chr2_list.txt")
else:
    generate_images(input_file, size=size, int_matrix=False, output_folder=folder_path, affix=("_" + chr_num), show_plot=False, no_norm=True, sequence_names_file="./chr2_list.txt")
