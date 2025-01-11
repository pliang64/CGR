import copy
import os
import sys
import time
import torch
import torch.nn as nn
import torch.utils.data.dataloader as dataloader
import torchvision
import torchvision.transforms as transforms
from torchvision.datasets import ImageFolder
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
import numpy as np
from PIL import Image
import random as r
from sklearn.metrics import confusion_matrix, balanced_accuracy_score, accuracy_score
from torch.autograd import Variable
import argparse
from cgr_by_chrom_v2 import generate_images
from collections import Counter



def parse_arguments():
    parser = argparse.ArgumentParser(description='Finetune and evaluate a model based on CGR plots. Run with no args for interactive commandline.')
    # Configuration parameter args
    parser.add_argument('--dataset_folder', type=str, default=dataset_folder, help='Dataset folder path')
    parser.add_argument('--save_folder', type=str, default=save_folder, help='Save folder path')
    parser.add_argument('--model_name', type=str, default=model_name, help='Model name')
    parser.add_argument('--weighted_labels', type=str, default=weighted_labels, help='Weighted labels')
    parser.add_argument('--seed', type=int, default=seed, help='Random seed')
    parser.add_argument('--chromnum', type=str, default=chromnum, help='Chr#. -1 for all chromosomes, 0 if whole genome.')
    # Hyper parameter args
    parser.add_argument('--batch_size', type=int, default=batch_size, help='Batch size')
    parser.add_argument('--shuffle', action='store_true', default=shuffle, help='Shuffle data')
    parser.add_argument('--epochs', type=int, default=epochs, help='Number of epochs')
    parser.add_argument('--lr', type=float, default=lr, help='Learning rate')
    parser.add_argument('--early_integration', action='store_true', help='Early Integration of all 19 chromosomes. Requires chromnum == -1.')
    parser.add_argument('--evaluate_ensemble', action='store_true', default=evaluate_ensemble, help='Evalaute ensemble of models. Requires a folder of models from saved_model.')
    # Classify mode
    parser.add_argument(
        '--classify',
        action='store_true',
        help='Classify an image using a single model, or folder of models for late integration voting. '
             'Makes use of dataset_folder for labels, model_name, saved_model, and cgr_target. '
             'When processing .fa files, save_folder is used to specify cgr output location.'
    )
    parser.add_argument('--saved_model', type=str, default=saved_model, help='Saved model path. If a folder is specified, late integration voting will be used with each of the models in the path.')
    parser.add_argument('--classify_target', type=str, default=classify_target, help='Classify folder or specific file (cgr/fa).') 
    
    
    return parser.parse_args()


def calc_class_weights():
    counts = [0] * num_classes
    tot = 0
    
    for _, labels in dataloaders["train"]:
        for label in labels:
            counts[label] += 1
            tot += 1
            
    weights = [(1 - (count / tot)) / (num_classes - 1) for count in counts]
    
    return torch.cuda.FloatTensor(weights)


def set_seed(x):
    r.seed(x)
    np.random.seed(x)
    torch.manual_seed(x)
    if torch.cuda.is_available(): 
        torch.cuda.manual_seed_all(x)
    torch.backends.cudnn.deterministic = True
    

#model training
def train_model(model,criterion,optimizer,scheduler,num_epochs,dataloader,dataset):
    since = time.time()

    best_model_wts = copy.deepcopy(model.state_dict())
    best_acc = 0.0

    for epoch in range(num_epochs):
        # Print conditions to reduce output
        if (epoch == 0 or 
                epoch == num_epochs - 1 or 
                (num_epochs < 100 and epoch % 5 == 0) or
                (num_epochs > 100 and epoch % 10 == 0) or
                num_epochs < 10
                # or True # Uncomment to always print
               ):
            print_output = True
        else:
            print_output = False
        if print_output:
            print("-"*10)
            print(f'Epoch {epoch}/{num_epochs-1}')
            print("-"*10)

        for phase in ['train','val']:
            if phase == 'train':
                model.train()
            else:
                model.eval()

            running_loss = 0.0
            running_corrects = 0

            
            y_preds = []
            y_true = []
            for inputs,labels in dataloader[phase]:
                inputs = inputs.to(device)
                labels = labels.to(device)
                outputs = model(inputs)
                _, preds = torch.max(outputs, 1)
                y_preds.extend(preds.cpu().numpy())
                y_true.extend(labels.cpu().numpy())
                optimizer.zero_grad()

                with torch.set_grad_enabled(phase == 'train'):
                    outputs = model(inputs)
                    _,preds = torch.max(outputs,1)
                    loss = criterion(outputs,labels)
                    if phase == 'train':
                        loss.backward()
                        optimizer.step()

                running_loss += loss.item() * inputs.size(0)
                running_corrects += torch.sum(preds == labels.data)
            if phase == 'train':
                scheduler.step()

            epoch_loss = running_loss / len(dataset[phase])
            epoch_acc = running_corrects.double() / len(dataset[phase])
            balanced_acc = balanced_accuracy_score(y_true, y_preds)
            
            
            if print_output:
                print(f'{phase} Loss: {epoch_loss:.4f} Acc: {epoch_acc:.4f}')
                print(f'{phase} Balanced acc: {balanced_acc:.4f}')

            # if phase == 'val' and epoch_acc > best_acc:
            if phase == 'val' and balanced_acc > best_acc:
                best_acc = balanced_acc
                best_model_wts = copy.deepcopy(model.state_dict())

        if print_output:
            print()

    time_elapsed = time.time() - since

    print(f'Training complete in {time_elapsed // 60:.0f}m {time_elapsed % 60:.0f}s')
    print(f'Best val acc: {best_acc:4f}')

    model.load_state_dict(best_model_wts)
    return model, best_acc

def eval_model(model, filename):
    model.to('cuda:0')
    y_preds = []
    y_true = []
    for inputs, labels in dataloaders["test"]:
        inputs = inputs.to('cuda:0')  # Move inputs to GPU
        labels = labels.to('cuda:0')  # Move labels to GPU

        outputs = model(inputs)
        _, preds = torch.max(outputs, 1)
        y_preds.extend(preds.cpu().numpy())  # Move predictions to CPU
        y_true.extend(labels.cpu().numpy())  # Move labels to CPU

    classes = train_dataset.classes

    cm = confusion_matrix(y_true, y_preds)
    df_cm = pd.DataFrame(cm, index=[i for i in classes], columns=[i for i in classes])

    balanced_acc = balanced_accuracy_score(y_true, y_preds)
    acc = accuracy_score(y_true, y_preds)  # Added this line for accuracy

    print("Balanced Accuracy:", balanced_acc)
    print("Accuracy:", acc)  # Added this line to print the accuracy

    plt.figure(figsize=(11, 10))
    plt.suptitle("Model Evaluation\nBalanced Accuracy: {:.4f}\nAccuracy: {:.4f}".format(balanced_acc, acc), fontsize=16, ha='center')
    plt.subplots_adjust(top=0.80)
    # plt.subplot(121)
    plt.title("Confusion Matrix")
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    sb.heatmap(df_cm, annot=True)
    
    # plt_path = os.path.join(save_folder,save_name+str_acc+".png")
    plt_path = filename+f"_Test{balanced_acc:.4f}.png"
    plt.tight_layout()
    plt.savefig(plt_path)
    print(f"Confusion matrix saved as {plt_path}")

    

# Interactive commandline functions
def prompt_update_var(message, default_value):
    user_input = input(f"{message} (default: {default_value}): ")
    return user_input if user_input else default_value

def update_config(config, var_list):
    # Print the menu of variables
    print("*"*10)
    print(f"Select a {config} variable to update:")
    for index, var_name in enumerate(var_list.keys(), start=1):
        print(f"{index}. {var_name}")
        
    while True:
        selection = input("Enter the number of the variable to update (or press Enter to finish): ")
        
        if not selection:
            break
        
        try:
            index = int(selection)
            if 1 <= index <= len(var_list):
                var_name = list(var_list.keys())[index - 1]
                new_value = prompt_update_var(f"Enter new {var_name}", var_list[var_name])
                var_list[var_name] = new_value
                print(f"Updated {var_name}: {var_list[var_name]}")
            else:
                print("Invalid selection. Please enter a valid number or press Enter to continue.")
        except ValueError:
            print("Invalid input. Please enter a valid number or press Enter to continue")
    return var_list


# Update datasets and dataloader using chr i
def update_dataset(i):
    global train_dataset, val_dataset, test_dataset, train_loader, val_loader, test_loader
    
    if i == 0:
        train_dataset = ImageFolder(root=data_path+"/training",transform=transform)
        val_dataset = ImageFolder(root=data_path+"/validation",transform=transform)
        test_dataset = ImageFolder(root=data_path+"/testing",transform=transform)
    else:
        train_dataset = ImageFolder(root=data_path+"/training/chr"+str(i),transform=transform)
        val_dataset = ImageFolder(root=data_path+"/validation/chr"+str(i),transform=transform)
        test_dataset = ImageFolder(root=data_path+"/testing/chr"+str(i),transform=transform)
    
    train_loader = dataloader.DataLoader(train_dataset,batch_size=batch_size,shuffle=shuffle)
    val_loader = dataloader.DataLoader(val_dataset,batch_size=len(val_dataset),shuffle=shuffle)
    test_loader = dataloader.DataLoader(test_dataset,batch_size=len(test_dataset),shuffle=shuffle)
    
    
def run_training_and_evaluation():
    # Train model
    global optimizer, str_acc
    train_size = len(train_dataset)
    model_current, acc = train_model(model,criterion,optimizer,torch.optim.lr_scheduler.OneCycleLR(optimizer,max_lr=lr,steps_per_epoch=round(train_size/batch_size), epochs=1,pct_start=0.99), 1,dataloaders,datasets)
    optimizer = torch.optim.Adam(model_current.parameters(), lr=lr/2)
    model_current, acc = train_model(model_current,criterion,optimizer,torch.optim.lr_scheduler.OneCycleLR(optimizer,max_lr=lr/2,steps_per_epoch=round(train_size/batch_size), epochs=epochs,pct_start=0.3), epochs,dataloaders,datasets)                          
    str_acc = f"{acc:.4f}"
    saved_model = os.path.join(save_folder,save_name+str_acc+".pkl")
    torch.save(model_current.state_dict(),saved_model)
    
    # Evaluate model
    model_current.load_state_dict(torch.load(saved_model))
    eval_model(model_current,saved_model)
    
def classify(model, image_path):    
    # Load and preprocess the image
    image = Image.open(image_path).convert("RGB")
    image = transform(image).unsqueeze(0).to(device)  # Add batch dimension
    
    # Make a prediction
    with torch.no_grad():
        outputs = model(image)
    
    # Process the prediction and print the result (modify as needed)
    _, predicted_class = torch.max(outputs, 1)
    # print("Predicted class:", class_labels[predicted_class.item()])
    return class_labels[predicted_class.item()]
    
def get_filelist(target, extension):
    filelist = []

    if os.path.isfile(target) and target.lower().endswith(extension):
        filelist.append(target)
    elif os.path.isdir(target):
        for file in os.listdir(target):
            if file.lower().endswith(extension):
                filelist.append(os.path.join(target, file))

    return filelist

def calculate_votes(votes_dict):
    # Find the class label with the highest vote count
    winning_class = max(votes_dict, key=votes_dict.get)
    
    # Calculate the percentage of votes for the winning class
    total_votes = sum(votes_dict.values())
    winning_percentage = (votes_dict[winning_class] / total_votes) * 100
    
    return winning_class, winning_percentage


def eval_ensemble(modellist):
    all_preds = []
    all_labels = []
    
    i=0;
    
    for model_path in modellist:
        i += 1
        print(f"{i}/{len(modellist)}: {model_path}")
        model.load_state_dict(torch.load(model_path))
        model.to(device)
        
        y_preds = []
        y_true = []
        
        for inputs, labels in dataloaders["test"]:
            inputs = inputs.to('cuda:0')
            labels = labels.to('cuda:0')
            
            with torch.no_grad():
                outputs = model(inputs)
                _, preds = torch.max(outputs, 1)
            
            y_preds.extend(preds.cpu().numpy())
            y_true.extend(labels.cpu().numpy())
        
        all_preds.append(y_preds)
        
        if len(all_labels) == 0:
            all_labels.extend(y_true)
    
    # Late Integration
    ensemble_preds = []
    for i in range(len(all_labels)):
        ith_preds = [y_preds[i] for y_preds in all_preds]
        most_common, _ = Counter(ith_preds).most_common(1)[0]
        ensemble_preds.append(most_common)
    
    classes = train_dataset.classes
    
    print(f"Test: {classes}")
    
    # Evaluate ensemble model
    cm = confusion_matrix(all_labels, ensemble_preds)
    df_cm = pd.DataFrame(cm, index=[i for i in classes], columns=[i for i in classes])
    
    balanced_acc = balanced_accuracy_score(all_labels, ensemble_preds)
    acc = accuracy_score(all_labels, ensemble_preds)
    
    print("Ensemble Balanced Accuracy:", balanced_acc)
    print("Ensemble Accuracy:", acc)
    
    # Confusion matrix
    plt.figure(figsize=(11, 10))
    plt.suptitle("Ensemble Model Evaluation\nBalanced Accuracy: {:.4f}\nAccuracy: {:.4f}".format(balanced_acc, acc), fontsize=16, ha='center')
    plt.subplots_adjust(top=0.80)
    plt.title("Confusion Matrix")
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    sb.heatmap(df_cm, annot=True)
    
    plt_path = f"ensemble_Test{balanced_acc:.4f}.png"
    plt.tight_layout()
    plt.savefig(plt_path)
    print(f"Confusion matrix saved as {plt_path}")
    
### Check device
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
torch.cuda.set_device(device)
torch.cuda.empty_cache()


    
valid_model_names = ["resnet18", "resnet101", "densenet121", "vit_b_16", "vgg11_bn"]

# Configuration parameters
dataset_folder = "../training_data_2"
save_folder = "../models"
model_name = "resnet18"
weighted_labels = "weighted"
seed=42
chromnum = "1"
early_integration = False
late_integration = False
classification = False
batch_classification = False
unprocessed_fa = False
evaluate_ensemble = False


saved_model = "../models/resnet18"
# saved_model = "../models/resnet18/chr2_resnet18_weighted_training_data_2_B16_E200_L0.01_BA0.4430.pkl"
# classify_target = "../input_images/SRR5891871_224_robust_gray_cgr.png"
classify_target = "../input_images/"
cgr_output_folder = classify_target

# Hyper parameters
batch_size = 16
shuffle = False
epochs = 200
lr = 0.01

# Transform for grayscale images
transform = transforms.Compose([
    transforms.Resize((224, 224)),
    transforms.Grayscale(num_output_channels=3),
    transforms.ToTensor(),
    transforms.Normalize(
        mean=[0.5, 0.5, 0.5],  # Use three-channel mean
        std=[0.5, 0.5, 0.5]   # Use three-channel std
    )
])

transform_vit = transforms.Compose([
    transforms.Resize((384, 384)),  # Make sure this matches your specific ViT variant
    transforms.Grayscale(num_output_channels=3),
    transforms.ToTensor(),
    transforms.Normalize(
        mean=[0.485, 0.456, 0.406],  # ImageNet mean
        std=[0.229, 0.224, 0.225]   # ImageNet std
    )
])

###################################################################################################

# Parse args if given, otherwise use an interactive script

if len(sys.argv) <2 or sys.argv[1] == "-f":
    # Interactive script
    print("Are you finetuning a model or classifying?\n"
                        "1.'finetune' (default)\n"
                        "2.'classify'\n"
                        "3.'ensemble'\n")
    while True:
        user_choice = input("Enter selection:")

        print(user_choice)

        # Finetune interactive script
        if user_choice.lower() == 'finetune' or user_choice == '' or user_choice == '1':
            # Finetune model
            print("Finetuning mode")
            print("Available models: "+', '.join(valid_model_names))
            
            config_vars = {
                "dataset_folder": dataset_folder,
                "save_folder": save_folder,
                "model_name": model_name,
                "weighted_labels": weighted_labels,
                "seed": seed,
                "chromnum": chromnum,
                "early_integration": early_integration
                }
            hyper_vars = {
                "batch_size": batch_size,
                "shuffle": shuffle,
                "epochs": epochs,
                "lr": lr
                }
            config_vars = update_config("Configuration Parameter", config_vars)
            hyper_vars = update_config("Hyper Parameter", hyper_vars)


            # Assign values from config_vars
            dataset_folder = config_vars["dataset_folder"]
            save_folder = config_vars["save_folder"]
            model_name = config_vars["model_name"]
            weighted_labels = config_vars["weighted_labels"]
            seed = int(config_vars["seed"])
            chromnum = config_vars["chromnum"]
            early_integration = bool(config_vars["early_integration"])

            # Assign values from hyper_vars
            batch_size = int(hyper_vars["batch_size"])
            shuffle = bool(hyper_vars["shuffle"])
            epochs = int(hyper_vars["epochs"])
            lr = float(hyper_vars["lr"])

            # Print the final vars before starting
            print("*"*10)
            print("Starting with the following settings:")
            for var_name, var_value in config_vars.items():
                print(f"{var_name}: {var_value}")
            for var_name, var_value in hyper_vars.items():
                print(f"{var_name}: {var_value}")
            print("*"*10)
            break  # Exit the loop after valid choice

        # Classify interactive script
        elif user_choice.lower() == 'classify' or user_choice == '2':
            # Call function for image classification
            print("Classification mode")
            classification = True
            class_vars = {
                "dataset_folder": dataset_folder,
                "model_name": model_name,
                "saved_model": saved_model,
                "seed": seed,
                "classify_target": classify_target
            }
            class_vars = update_config("Classification Parameter", class_vars)

            dataset_folder = class_vars["dataset_folder"]
            model_name = class_vars["model_name"]
            saved_model = class_vars["saved_model"]
            seed = class_vars["seed"]
            classify_target = class_vars["classify_target"]
            if os.path.isdir(classify_target):
                late_integration = True
                class_vars["late_integration"] = late_integration
            if os.path.isdir(classify_target):
                batch_classification = True
                class_vars["batch_classification"] = batch_classification
            folder_contains_fa = os.path.isdir(classify_target) and any(file.lower().endswith(".fa") for file in os.listdir(classify_target))
            if folder_contains_fa or classify_target.lower().endswith(".fa"):
                unprocessed_fa = True
                cgr_output_folder = prompt_update_var(f"Specify a folder for output cgr files:", cgr_output_folder)
                class_vars["cgr_output_folder"] = cgr_output_folder
            
            
            print("*"*10)
            print("Starting with the following settings:")
            for var_name, var_value in class_vars.items():
                print(f"{var_name}: {var_value}")
            print("*"*10)
            break  # Exit the loop after valid choice
        
        # Evaluate ensemble interactive script
        elif user_choice.lower() == 'ensemble' or user_choice == '3':
            # Call function for image classification
            print("Ensemble Evaluation mode")
            evaluate_ensemble = True
            ens_vars = {
                "dataset_folder": dataset_folder,
                "model_name": model_name,
                "saved_model": saved_model
            }
            ens_vars = update_config("Ensemble Parameter", ens_vars)

            dataset_folder = ens_vars["dataset_folder"]
            model_name = ens_vars["model_name"]
            saved_model = ens_vars["saved_model"]
            break  # Exit the loop after valid choice
            
        else:
            print("Invalid choice. Please enter either 'finetune' or 'classify'.")


# Parse args    
else:
    args = parse_arguments()
    
    if args.classify:
        dataset_folder = args.dataset_folder
        model_name = args.model_name
        classify_target = args.classify_target
        saved_model = args.saved_model
        cgr_output_folder = args.save_folder
    elif args.evalaute_ensemble:
        dataset_folder = args.dataset_folder
        model_name = args.model_name
        saved_model = args.saved_model  
    else:
        dataset_folder = args.dataset_folder
        save_folder = args.save_folder
        model_name = args.model_name
        weighted_labels = args.weighted_labels
        seed = int(args.seed)
        chromnum = args.chromnum
        batch_size = int(args.batch_size)
        shuffle = bool(args.shuffle)
        epochs = int(args.epochs)
        lr = float(args.lr)
    


    
if weighted_labels.lower() == 'weighted':
    weighted_ce = True
else:
    weighted_labels = "unweighted"
    weighted_ce = False


set_seed(seed)

###################################################################################################
# Prepare Dataloaders
global train_dataset, val_dataset, test_dataset, train_loader, val_loader, test_loader
data_path = dataset_folder

if int(chromnum) < 0 and early_integration:
    train_subdirs = [os.path.join(data_path,"training",chrm) for chrm in os.listdir(os.path.join(data_path,"training")) if os.path.isdir(os.path.join(data_path,"training",chrm))]
    val_subdirs = [os.path.join(data_path,"val",chrm) for chrm in os.listdir(os.path.join(data_path,"val")) if os.path.isdir(os.path.join(data_path,"val",chrm))]
    test_subdirs = [os.path.join(data_path,"test",chrm) for chrm in os.listdir(os.path.join(data_path,"test")) if os.path.isdir(os.path.join(data_path,"test",chrm))]
    
    train_dataset = []
    val_dataset = []
    test_dataset = []
    
    if model_name == "vit_b_16":
        transform = transform_vit
    
    # Load all chromosomes for early integration
    for train_subdir in train_subdirs:
        tmp_dataset = ImageFolder(root=train_subdir, transform=transform)
        train_dataset.append(tmp_dataset)
    for val_subdir in val_subdirs:
        tmp_dataset = ImageFolder(root=val_subdir, transform=transform)
        val_dataset.append(tmp_dataset)
    for test_subdir in test_subdirs:
        tmp_dataset = ImageFolder(root=test_subdir, transform=transform)
        test_dataset.append(tmp_dataset)
        
    train_loader = dataloader.DataLoader(train_dataset,batch_size=batch_size,shuffle=shuffle)
    val_loader = dataloader.DataLoader(val_dataset,batch_size=len(val_dataset),shuffle=shuffle)
    test_loader = dataloader.DataLoader(test_dataset,batch_size=len(test_dataset),shuffle=shuffle)
elif int(chromnum) == 0:
    update_dataset(0)
else:
    if int(chromnum) < 0:
        update_dataset(1)
    else:
        update_dataset(chromnum)

datasets = {"train":train_dataset, "val":val_dataset, "test":test_dataset}
dataloaders = {"train":train_loader,"val":val_loader, "test":test_loader}

num_classes = len(train_dataset.classes)
class_labels = train_dataset.classes

# Calculate the total number of images
total_images = len(train_dataset)
if not classification:
    print(f"Total number of classes: {num_classes}")
    print(f"Total number of images: {total_images}")

# Set up models
if model_name == "resnet18":
    resnet_model = "resnet18-5c106cde"
    model_path = "../resnet/"+resnet_model+".pth"     # Path to the locally saved ResNet model weights
    save_folder = os.path.join(save_folder, model_name+"_earlyIntegration" if early_integration else model_name)
    os.makedirs(save_folder, exist_ok=True)
    model = torchvision.models.resnet18(weights=None,progress=True)
    model.load_state_dict(torch.load(model_path))
    model.to(device)
    model.fc = nn.Linear(in_features=512,out_features=num_classes,bias=True).to(device)
elif model_name == "resnet101":
    resnet_model = "resnet101-5d3b4d8f"
    model_path = "../resnet/"+resnet_model+".pth"     # Path to the locally saved ResNet model weights
    save_folder = os.path.join(save_folder, model_name+"_earlyIntegration" if early_integration else model_name)
    os.makedirs(save_folder, exist_ok=True)
    model = torchvision.models.resnet101(weights=None,progress=True)
    model.load_state_dict(torch.load(model_path))
    model.to(device)
    model.fc = nn.Linear(in_features=2048,out_features=num_classes,bias=True).to(device)
elif model_name == "densenet121":
    densenet_model = "densenet121"
    model_path = "../densenet/"+densenet_model+".pth"     # Path to the locally saved Densenet model weights
    save_folder = os.path.join(save_folder, model_name+"_earlyIntegration" if early_integration else model_name)
    os.makedirs(save_folder, exist_ok=True)
    model = torchvision.models.densenet121(weights=None,progress=True)
    model.load_state_dict(torch.load(model_path))
    model.to(device)
    model.fc = nn.Linear(in_features=2048,out_features=num_classes,bias=True).to(device)
elif model_name == "vit_b_16":
    # print(model)
    vit_model = "vit_b_16-c867db91"
    model_path = "../vit/"+vit_model+".pth"     # Path to the locally saved ViT model weights
    save_folder = os.path.join(save_folder, model_name+"_earlyIntegration" if early_integration else model_name)
    os.makedirs(save_folder, exist_ok=True)
    model = torchvision.models.vit_b_16(weights=None,progress=True)
    model.load_state_dict(torch.load(model_path))
    model.to(device)
    model.heads = nn.Sequential(nn.Linear(768,19,bias=True)).to(device)
elif model_name == "vgg11_bn":
    # print(model)
    vgg_model = "vgg11_bn-6002323d"
    model_path = "../vgg/"+vgg_model+".pth"     # Path to the locally saved vgg model weights
    save_folder = os.path.join(save_folder, model_name+"_earlyIntegration" if early_integration else model_name)
    os.makedirs(save_folder, exist_ok=True)
    model = torchvision.models.vgg11_bn(weights=None,progress=True)
    model.load_state_dict(torch.load(model_path))
    model.to(device)
    model.classifier[6] = nn.Linear(4096,19,bias=True).to(device)
else:
    raise ValueError(f"Invalid model name '{model_name}'. Valid options are: {', '.join(valid_model_names)}")

# Define Optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=lr)
if weighted_ce:
    criterion = nn.CrossEntropyLoss(weight=calc_class_weights()).to(device)
else:
    criterion = nn.CrossEntropyLoss().to(device)

chr_string = "" if int(chrnum) == 0 else f"chr{chrnum}_"
# save_name = "chr"+ chromnum + "_" + model_name+"_"+weighted_labels+"_"+os.path.basename(data_path)+"_"
save_name = f"{chr_string}{model_name}_{weighted_labels.lower()}_{os.path.basename(data_path)}_B{batch_size}_E{epochs}_L{lr}_BA"


str_acc = None

if __name__ == '__main__':
    # # Evaluate model
    # saved_model = "../models/resnet18/chr1_resnet18_weighted_training_data_2_B16_E200_L0.01_BA0.4079"
    # model.load_state_dict(torch.load(saved_model+".pkl"))
    # eval_model(model,saved_model)
    
    if classification and evaluate_ensemble:
        print("Error: Select either classification, finetuning, or evaluate_ensemble.")
        sys.exit(1)  # The argument (1) is a status code indicating an error occurred
    
    # Classification
    if classification:
        filelist = []
        # Get png files in classify_target
        filelist = get_filelist(classify_target,".png")
        # Process fa into CGR if necessary
        if unprocessed_fa:
            falist = get_filelist(classify_target,".fa")
            for sequence_file in falist:
                filelist += generate_images(sequence_file, size=224, output_folder=cgr_output_folder, robust_scaler=True)
        # print(filelist)
        
        # Classify
        modellist = get_filelist(saved_model,".pkl")
        # print(modellist)
        
        loaded_models = []
        
            
        # print(len(loaded_models))
        # print(loaded_models)
        
        for image in filelist:
            #dict
            chr_predictions = {}
            
            print(f"Classifying: {image}")
            
            for load_model in modellist:
                model.load_state_dict(torch.load(load_model))
                model.eval().to(device)
                # print(load_model)
                predicted_class = classify(model, image)
                if predicted_class in chr_predictions:
                    chr_predictions[predicted_class] += 1
                else:
                    chr_predictions[predicted_class] = 1
                
            # print(chr_predictions)
            winning_class, winning_percentage = calculate_votes(chr_predictions)
            if late_integration:
                print(f"Vote Winner: {winning_class}, {winning_percentage}%")
            else:
                print(f"Predicted Class: {winning_class}") 
    
    # Evaluate ensemble
    elif evaluate_ensemble:
        modellist = get_filelist(saved_model,".pkl")
        print(saved_model)
        print(f"Models: {modellist}")
        # print(f"length: {len(modellist)}")
        eval_ensemble(modellist)
    
    # Finetune
    else:
        if int(chromnum) < 0 and not early_integration:
            training_subdirs = sorted([name for name in os.listdir(os.path.join(data_path,"training")) if os.path.isdir(os.path.join(os.path.join(data_path,"training"), name))], key=lambda name: int(name[3:]))
            print(training_subdirs)
            print(f"Training 1-{len(training_subdirs)}")
            for i in range(1, len(training_subdirs)+1):
            # for i in range(1, 20): # 1-19 inclusive
                print("*"*20)
                print(f"Processing Chromosome {i}/{len(training_subdirs)}: {data_path}/training/chr{i}")
                save_name = f"chr{i}_{model_name}_{weighted_labels.lower()}_{os.path.basename(data_path)}_B{batch_size}_E{epochs}_L{lr}_BA"
                update_dataset(i)
                run_training_and_evaluation()
        else:
            run_training_and_evaluation()

   
                # pass
            ################################
                
#         model.load_state_dict(torch.load(saved_model))
#         model.eval()
#         device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
#         model.to(device)

#         for file in filelist:
#             classify(model, file)

# Square the confusion matrix
# Add accuracy (along with balance) and test accuracy to confusion matrix
# Try 224 no-normalization, 512 resolution, 512 resolution with no-normalization


# ../models_test
# /home/avvu/scratch/training_data_1024