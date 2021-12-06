# MERIDA Snakemake pipeline

Steps to running MERIDA on CCLE dataset 
1. Install packages needed for R script Make_Bimodality_Input.R
2. run all lines of the script. this will produce two files: AAC_file.txt and Input_Matrix.txt. Download these files to your machine 
3. I ran MERIDA on the Cedar cluster of Compute Canada. The steps to do so are the following

    a. Install CPLEX 12.8.0 into your compute canada account. I found the steps at this link (https://docs.computecanada.ca/wiki/CPLEX/en) easy to follow 
    
    b. clone this repository into your scratch directory. Also, copy your AAC file and Input Matrix to Compute Canada
  
    e. Unzip the MERIDA.tar.gz file to access the MERIDA directory
  
    f. navigate into the MERIDA directory 
  
    g. In the file `Example_Config1.txt`, replace the paths to the Input Matrix and the AAC file to be the paths to those files in your account.
  
    h. Run the command `./build/MERIDA_ILP Example_Config1.txt no` to train the model. 
  
    i. For usage instructions, run `./build/MERIDA_ILP` with no parameters.
  
    j. If this command does not work and compilation is required, do the following:
  
    k. navigate into the build directory and remove its contents 
      - run `cmake` (ensure the results are stored in the build directory) 
      - navigate to the parent directory and run `make` to generate the exectuable
