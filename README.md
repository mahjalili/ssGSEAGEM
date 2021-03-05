# ssGSEAGEM

<h3>Introduction</h3>
The ssGSEAGEM introduced a new conceptual framework for transcriptional data integration into a genome-scale metabolic model, to accurately predict the metabolic phenotypes of a model. (such as <i>Saccharomyces cerevisiae</i>).<br>
This framework could accurately predict metabolic phenotypes and outperform the accuracy. We believe that this framework could be adopted almost to all data integration algorithms and be enhanced their prediction power.<br>
![ssGSEAGEM workflow](https://github.com/mahjalili/ssGSEAGEM/blob/main/Images/workflow.jpg?raw=true)
<br>
<br>
<h3>Usage</h3>
<b>Step 1: </b>First prepare conventional gene expression data and convert gene expression levels to reaction levels using GPR associations (Script: Matlab/datasets/GIMME).<br>
<b>Step 2: </b>Using first step data and geneset of GEM model (Script: Matlab/ssGSEA/geneset) prepar ssGSEA scores of genesets (Script: R script).<br>
<b>Step 3: </b>Using first step data and gssGSEA scores prepare ssGSEA-GIMME scores of genesets and convert gene expression levels to reaction levels using GPR associations (Script: Matlab/datasets/ssGSEA).<br>
b>Step 4: </b>Make model by integrate GEM model (Matlab/models) and reactions expression levels provided by previous steps (Script: Matlab/makemodel).
<br>
<br>
<h3>Citation</h3>
Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Bas Teusink, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi, A conceptual framework for transcriptional data integration into a genome-scale metabolic model. 2021.
