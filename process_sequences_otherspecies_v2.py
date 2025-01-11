from cgr_by_chrom_v3 import generate_images
from process_folders_species import process_files
import os
import sys
import fnmatch
from matplotlib.colors import LogNorm

if len(sys.argv) != 5:
    print(f"Usage: python {sys.argv[0]} <chr#> <size N> <no_norm/robust_scaler> <Vitis_vinifera T>")
    sys.exit(1)  # Exit
    
def search_files(input_folder, pattern):
    matching_files = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if sys.argv[4] == "True":
                if fnmatch.fnmatch(file, pattern) and "Vitis_vinifera" in file:
                    file_path = os.path.join(root, file)
                    matching_files.append(file_path)
            elif sys.argv[4] == "False":
                if fnmatch.fnmatch(file, pattern) and "Vitis_vinifera" not in file:
                    file_path = os.path.join(root, file)
                    matching_files.append(file_path)
            else:
                print(f"Error: argv[4] must be True or False")
                sys.exit(1)  # Exit
                
    return matching_files

input_folder = "/home/avvu/scratch/data/species"
# pattern = "*."+"2"+".fa"  # Example pattern with wildcard
pattern = "*."+sys.argv[1]+".fa"  # Example pattern with wildcard
# out_folder = "./test"
out_folder = f"/home/avvu/scratch/cgr/species_{sys.argv[2]}_{sys.argv[3]}"
os.makedirs(out_folder, exist_ok=True)


# test_input_files = [os.path.join(input_folder,"Vitis_chungii.2.fa")]


# Generate folders for each of the files in input_folder
process_files(input_folder, out_folder)

# Process each file in input_folder
input_files = search_files(input_folder, pattern)
# for sequence_file in test_input_files:
for sequence_file in input_files:

    # Check by chrom
    filename = os.path.basename(sequence_file)
    prefix = filename.split(".")[0]
    #prefix = prefix.split("_")[0]
    chr_num = filename.split(".")[1]
    chr_num = "chr" + chr_num
    
    # Set folder path
    folder_path = os.path.join(out_folder, prefix, chr_num, "")
    print(folder_path)
    #print(folder_path)
    
    # Generate folder for the chrom of the sample if it does not exist
    if not os.path.exists(folder_path):
        os.makedirs(folder_path, exist_ok=True)
        print(f"Folder '{folder_path}' created.")
        
    # CGR
    print(f"Processing '{filename}':")
    print(folder_path)
    if sys.argv[3] == "no_norm":
        generate_images(sequence_file, size=int(sys.argv[2]), int_matrix = False, output_folder=folder_path, affix=("_" + chr_num), show_plot=True, no_norm=True)
    elif sys.argv[3] == "robust_scaler":
        generate_images(sequence_file, size=int(sys.argv[2]), int_matrix = False, output_folder=folder_path, affix=("_" + chr_num), show_plot=True, robust_scaler=True)
