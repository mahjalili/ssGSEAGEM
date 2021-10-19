# ssGSEAGEM

<h3>Introduction</h3>
The ssGSEAGEM introduced a new conceptual framework for transcriptional data integration into a genome-scale metabolic model, to accurately predict the metabolic phenotypes of a model. (such as <i>Saccharomyces cerevisiae</i>).<br>
This framework could accurately predict metabolic phenotypes and outperform the accuracy. We believe that this framework could be adopted almost to all data integration algorithms and be enhanced their prediction power.<br><br>
<img src="https://github.com/mahjalili/ssGSEAGEM/blob/main/Images/workflow.png?raw=true" alt="ssGSEAGEM workflow">
<br>
<h3>Usage</h3>
<b>Step 1: </b>First prepare conventional gene expression data and convert gene expression levels to reaction levels using GPR associations (Script: Matlab/datasets/standard.m).<br>
<b>Step 2: </b>Using first step data and geneset of GEM model (Script: Matlab/ssGSEA), prepare ssGSEA scores of genesets (Script: R scripts).<br>
<b>Step 3: </b>Using first step data and ssGSEA scores, prepare ssGSEA-GIMME scores of genesets and convert gene expression levels to reaction levels using GPR associations (Script: Matlab/datasets/ssGSEA).<br>
<b>Step 4: </b>Make models by integrate GEM model (Matlab/models) and reaction expression levels provided by previous steps (Script: Matlab/makemodel).
<br>
<br>
<h3>Citation</h3>
<b>Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi, Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models. 2021.</b>
