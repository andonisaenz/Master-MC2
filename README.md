Master Thesis
---------------
This repository contains all the code implemented in my project. 

1. The file `Main_script.r` contains all the code implemented in the RNAseq analysis pipeline, including:
    * All libraries that were used in the project.
    * Reference genome indexing and reads alignment.
    * Generation of counts matrix.
    * Data pre-processing and filtering steps.
    * Differential Expression Analysis with DESeq2.
    * Data normalization methods and visualization.
    * Obtention of differentially expressed gene lists.
    * Functional analysis with ORA and GSEA.
      
2. The file `TCR_analysis` contains all the code implemented in the TCR clonality analysis pipeline, including:
    * Extraction of TCR clonal repertoire information with MiXCR (this part was executed in the Windows Powershell, not in R).
    * TCR clonotype analysis and visualization with _immunarch_.
    * Analysis of TCR clonal repertoire diversity.
   
3. The file `MultiQC_analysis` contains all the code implemented to generate plots for the outputs of the quality control analysis with _MultiQC_.
