from cgr_by_chrom_v3 import generate_images
from process_folders_species import process_files
import os
import sys
import fnmatch
from matplotlib.colors import LogNorm

def search_files(input_folder, pattern):
    matching_files = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if fnmatch.fnmatch(file, pattern) and "Vitis_vinifera" not in file:
                file_path = os.path.join(root, file)
                matching_files.append(file_path)
    return matching_files

input_folder = "/home/avvu/projects/rrg-pliang-ac/lianglab-rrg/grape/data4Andrew/4400samples/vv_cultivars/fasta"
# out_folder = "./test"
out_folder = f"/home/avvu/scratch/cgr/vvcultivar_{sys.argv[2]}_{sys.argv[3]}"

sequence_files = [os.path.join(input_folder,sys.argv[1])+".fa"]
# sequence_files = ["/home/avvu/projects/rrg-pliang-ac/lianglab-rrg/grape/data4Andrew/4400samples/Pinot_samplist.chr5.fa"]
print(f"Processing: {sys.argv[1]}")


# test_input_files = [os.path.join(input_folder,"erbamat.fa")]
# print(test_input_files)

# Generate folders for each of the files in input_folder
# process_files(input_folder, "/home/avvu/scratch/cgr/species_nonorm")

# Process each file in input_folder
# input_files = search_files(input_folder, pattern)
for sequence_file in sequence_files:
# for sequence_file in test_input_files:
    # Check by chrom
    filename = os.path.basename(sequence_file)
    prefix = filename.split(".")[0]
    #prefix = prefix.split("_")[0]
    # chr_num = filename.split(".")[1]
    # chr_num = "chr" + chr_num
    
    # Set folder path
    folder_path = os.path.join(out_folder, prefix,"")

    
    # Generate folder for the chrom of the sample if it does not exist
    if not os.path.exists(folder_path):
        os.makedirs(folder_path, exist_ok=True)
        print(f"Folder '{folder_path}' created.")
        
    # CGR
    print(f"Processing '{filename}':")
    print(folder_path)
    if sys.argv[3] == "no_norm":
        generate_images(sequence_file, size=int(sys.argv[2]), int_matrix = False, output_folder=folder_path, no_norm=True)
    elif sys.argv[3] == "robust_scaler":
        generate_images(sequence_file, size=int(sys.argv[2]), int_matrix = False, output_folder=folder_path, robust_scaler=True)
