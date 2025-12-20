This is data and code for the paper "Exact Coherent Structures of Sheared Double-Diffusive Convection".

Because the storage limit, we do not provide data for making supplementary videos. To avoid Windows path-length limit (Error 0x80010135) when extracting zip file, we highly recommend downloading data via Github clone.

Also, you need to install Git LFS (used to save large files on Github) to load full large files from github to your computer after cloning.
```bash
# step 1: clone data, large files will be saved by pointers instead of real files. [run this in WSL terminal, via VScode]
git clone https://github.com/duc-nguyen1410/Exact_Coherent_Structures_of_Sheared_Double_Diffusive_Convection.git

# step 2: get full data of large files [run this in Git bash, via VScode]
git lfs install # install if needed
cd Exact_Coherent_Structures_of_Sheared_Double_Diffusive_Convection
git lfs pull  # this downloads the actual large files
```

Running Jupyter notebooks, it may requires some libraries:
```bash
# create conda environment with libraries within anaconda
conda create -n test-env python=3.13 numpy matplotlib scipy pandas jupyter ipykernel pandas netCDF4

# or install libraries 
pip install numpy matplotlib scipy pandas jupyter ipykernel pandas netCDF4

## please select exact python kernel when running notebook. 
# for example, here select kernel "test-env", or use your available python kernel
```