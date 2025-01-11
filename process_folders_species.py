import os

def process_files(input_folder, out_folder=""):
    prefix_set = set()

    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith("fa"):
                # Extract the prefix
                prefix = file.split(".")[0]
                #prefix = prefix.split("_")[0]
                #print(file)

                # Add the prefix to the set (removes duplicates)
                prefix_set.add(prefix)

                # Rename the file
                #old_path = os.path.join(root, file)
                #new_path = os.path.join(root, prefix)
                #os.rename(old_path, new_path)
    
    # Remove entries starting with "merged"
    prefix_set = {item for item in prefix_set if not item.startswith("merged")}

    # Create folders in the current directory
    current_directory = os.getcwd()
    for prefix in prefix_set:
        folder_path = os.path.join(current_directory, out_folder, prefix)
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)
            print("Created folder:", folder_path)
#         else:
#             print("Folder already exists:", folder_path)
            
# Example usage
#input_folder = "/home/avvu/projects/rrg-pliang-ac/lianglab-rrg/grape/data4Andrew/4400samples"
#process_files(input_folder)
