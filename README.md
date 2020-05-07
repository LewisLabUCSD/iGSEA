# Bioinformatics_Tools
GSEA analysis on DESeq2 results (Future Plan: Name conversion, database query, ...)

## Quick Start 

0. *Download this package* on your machine.

1. You can run the examples now: 
   - Example 1. - *'Run_gsea_ex1.r'* 
     This example is designed to demonstrate how to run gsea analysis on the input DESeq2 result files (*'.csv'*) that are all put in the same directory. 

   - Example 2. - *'Run_gsea_ex2.r'* 
     This example is designed to demonstrate how to run gsea analysis on the input DESeq2 object: a *'.rda'* file  
 
2. To get start on your own dataset, you need to:
   - a) create a /subfolder inside the */Data* and copy your DESeq2 files to the /Data/subfolder folder.
   - b) setting the parameters in the section of *'Setting Parameter'* (see the details in the above two examples).
   - c) Note that, the GSEA has abundant human gene sets (and a few sets of mouse genes). Thus, we encourge users to convert their transcript IDs to *'human gene names'* (see below for the ID mapping procedures).
   
3. ID conversion: to sucessfully use this package, there are two ways:
   - a) users need to have their DESeq2 results contain a column named *"gene_name"* to indicate the human gene name for your transcript.
   - b) users can also provide an IDmapping file and assign it to the parameter *"IDmappingfile"* in your own script, by which our program will automatically convert your transcript ID to human gene name (please see the Tutorial *"Tutorial on obtaining IDmapping file.pdf"* for obtaining an ID-mapping file).
   
   

  
 
