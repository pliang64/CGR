import os
import shutil
import random
import argparse

def split_data(input_folder, output_folder, chrom, seed, log_path, num_file_limit):
    
    # Define the distribution for training, testing and validation data
    distribution_config = {
        3: (1, 1, 1),
        4: (2, 1, 1),
        5: (3, 1, 1),
        6: (4, 1, 1),
        7: (4, 1, 2),
        8: (5, 1, 2),
        9: (6, 1, 2)
    }
    train_dist = 0.7 # 0-train_dist will be distributed to training
    test_dist = 0.2 + train_dist # train_dist-test_dist will be distributed to testing
    # Remining 1-test_dist will be for validation
    
    chrom_folder = "" if chrom == 0 else f"chr{chrom}"

        
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    subfolders = [f.name for f in os.scandir(input_folder) if f.is_dir()]

    for species in subfolders:
        subfolder_path = os.path.join(input_folder, species, chrom_folder)
        if not os.path.isdir(subfolder_path):
            continue

        filenames = [f.name for f in os.scandir(subfolder_path) if f.is_file() and f.name.endswith('.png')]
        random.seed(seed)
        random.shuffle(filenames)

        num_files = len(filenames)
        max_num_files = min(num_file_limit, num_files)  # Process up to the first 60 files

        if max_num_files < 10:
            if max_num_files in distribution_config:
                num_train, num_val, num_test = distribution_config[max_num_files]
                train_files = filenames[:num_train]
                test_files = filenames[num_train:num_train + num_test]
                val_files = filenames[num_train + num_test:]
            else:
                with open(log_path, "a+") as file:
                    unsupported_message = f"{species}\n"
                    if unsupported_message not in file.read():
                        file.write(unsupported_message)
                    continue  # Skip the current iteration of the loop
        else: # if 10 < max_num_files < num_file_limit
            train_end = int(train_dist * max_num_files)
            test_end = int(test_dist * max_num_files)

            train_files = filenames[:train_end]
            test_files = filenames[train_end:test_end]
            val_files = filenames[test_end:num_file_limit]

        output_training = os.path.join(output_folder, 'training', chrom_folder, species)
        output_testing = os.path.join(output_folder, 'testing', chrom_folder, species)
        output_validation = os.path.join(output_folder, 'validation', chrom_folder, species)
        os.makedirs(output_training, exist_ok=True)
        os.makedirs(output_testing, exist_ok=True)
        os.makedirs(output_validation, exist_ok=True)
        
        
        for filename in train_files:
            src_path = os.path.join(subfolder_path, filename)
            dest_path = os.path.join(output_training, filename)
            try:
                os.symlink(os.path.abspath(src_path), dest_path)
            except FileExistsError as e:
                print(f"Symlink creation failed: {e}")

        for filename in test_files:
            src_path = os.path.join(subfolder_path, filename)
            dest_path = os.path.join(output_testing, filename)
            try:
                os.symlink(os.path.abspath(src_path), dest_path)
            except FileExistsError as e:
                print(f"Symlink creation failed: {e}")

        for filename in val_files:
            src_path = os.path.join(subfolder_path, filename)
            dest_path = os.path.join(output_validation, filename)
            # print(f"'{dest_path}' -> '{src_path}'")
            try:
                os.symlink(os.path.abspath(src_path), dest_path)
            except FileExistsError as e:
                print(f"Symlink creation failed: {e}")

def check_chrom_paths(folder_paths):
    existing_folders = [folder for folder in folder_paths if os.path.exists(folder)]
    
    if existing_folders:
        print("Warning: The following folders already exist:")
        for folder in existing_folders:
            print(f"  {folder}")
        
        user_input = input("Do you want to delete them? (Y/n): ").strip().lower()
        if user_input == 'y' or user_input == '':
            for folder in existing_folders:
                try:
                    shutil.rmtree(folder)  # Force delete the folder and its contents
                    print(f"'{folder}' folder deleted.")
                except OSError as e:
                    print(f"Error: Unable to delete '{folder}'. Reason: {e}")
                    return False
            return True
        else:
            print(f"Skipping current folders.")
            return False
    else:
        #print("No existing folders found.")
        return True
                

def main():
    parser = argparse.ArgumentParser(description='Split and organize labeled training data')
    parser.add_argument('--input', required=True, help='Path to input folder containing labeled training data')
    parser.add_argument('--output', required=True, help='Path to output folder for organized data')
    parser.add_argument('--chrom', type=int, default=-1, help='Specify chromosome (1-19), -1 to process every chromosome, or 0 for whole genomes.')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for data splitting')
    parser.add_argument('--max_files', type=int, default=60, help='Maximum number of files per label')
    args = parser.parse_args()
    
    if args.chrom == -1:
        loop_range = range(1, 20)
    elif 1 <= args.chrom <= 19:
        loop_range = range(args.chrom, args.chrom + 1)
    elif args.chrom == 0:
        print("Processing whole genome training data")
    else:
        print("Error: Invalid chromosome value. Please provide a value between 1 and 19 or -1.")
        return
    

        print(f"Performing loop for chromosome {i}")
    
    # Excluded log (fewer than 3 samples)
    log = "excluded_species.txt"
    log_path = os.path.join(args.output, log)
    if os.path.exists(log_path):
        os.remove(log_path)
        
        
    subfolders = ["testing", "training", "validation"]
    
    if args.chrom == 0:
        # Whole genomes
        folder_paths = [os.path.join(args.output, subfolder) for subfolder in subfolders]
        path_check = check_chrom_paths(folder_paths)
        if path_check:
            print(f"Processing.")
            split_data(args.input, args.output, 0, args.seed, log_path, args.max_files)
        else:
            print(f"Skipping.")
    else:    
        for i in loop_range:
            # Check if chr was already processed
            folder_paths = [os.path.join(args.output, subfolder, "chr"+str(i)) for subfolder in subfolders]
            path_check = check_chrom_paths(folder_paths)
            if path_check:
                print(f"Processing chr'{i}'")
                split_data(args.input, args.output, i, args.seed, log_path, args.max_files)
            else:
                print(f"Skipping chr'{i}'")
                continue

    
if __name__ == '__main__':
    main()
