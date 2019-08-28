# GPECOC

Genetic Programming based Error Correcting Output Codes

This is the implementation for paper: [A Novel Error-Correcting Output Code Algorithm Based on Genetic Programming](<https://www.sciencedirect.com/science/article/abs/pii/S2210650219302044>)

## Acknowledgement

- Codes about Genetic Programming is modified from  [Pyevolve](<https://github.com/perone/Pyevolve>)

- Codes about Output-Code-Classifier is modified from  [scikit-learn 0.18](<https://github.com/scikit-learn/scikit-learn/tree/0.18.X>)

## Environment

- **Windows 10 64 bit** 

- **python 2.7**

- **Excel**

  Enable the macro in excel, so it can extract result from file automatically.

- **scikit-learn 0.18**

  Anaconda is strongly recommended, run the following command in the Powershell, all necessary python packages for this project will be installed:

  ```shell
  conda install scikit-learn==0.18
  ```

## Dataset

- **Data format**

  Data should be put into the folder ```gpecoc/data```. Each dataset should be divided into "dataname_train.data", "datanam_test.data" and "dataname_validation.data". In the sub-dataset, each column is a sample, the first line represents the labels, the rest are feature space. There are two examples datasets in the folder ```gpecoc/data```. Please note that invalid sample, such as value missed, will cause errors.

- **Data processing**

  Feature Selection and Scaling will be done automatically. 


## Multi-processing Mode

- **Config**

  Firstly, make configuration in ```Configurations.py```. In Multi-processing mode, you need not to pay attention to 'dataName' and 'aimFolder'.
  
- **Run the following command**

  It will traversal all datasets given by the main function in ```_ParallelRunner.py```, each dataset will be run for 10 times.

  ```python
  python _ParallelRunner.py
  ```

- **Analyze result**

  All result infos will be written into the folder. For example, if you set version = "8.80", result infos will be found in 
  ($root_path)/Results/8.80/

## Single-processing Mode

This is useful when you want to debug.

- **Config**

  Make configuration in ```Configurations.py```. In Single-processing mode, 'dataName' and 'aimFolder' should also be set.
  
- **Run the following command**

  It will do training and testing on the dataset given in the ```Configurations.py```.

  ```python
  python AGpEcocStart.py
  ```

- **Parse and analyze result**

   In this Mode, part of the result will be printed on the terminal. You can find all result in the Results folder.
   But, there will be no automatic analyzing.

## Attention

Make sure the version folder is not exit every time you run it. Or errors will happen in the following way:
- In Single-processing Mode, the result in the terminal is right, but the result written ino the file might be wrong.
- In Multi-processing Mode, the result could not be read and parsed correctly.

A suggestion is to change the 'version' in the ```Configurations.py``` every time you run it.
