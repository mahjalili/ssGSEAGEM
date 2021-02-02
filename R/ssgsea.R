################################################################################################################
## Filename: ssgsea.r
## Created: 2020-08-26 19:14:31
## Author(s): Modified by Mahdi Jalili (origin by Karsten Krug)
##
## Purpose: 
##      - Wrapper to ssGSEA script to perform single sample Gene Set Enrichment analysis.
##          
## Instructions:  
##      - Source the script into a running R-session:
##          - RStudio: open file and press 'Source' in the upper right part of the editor window 
##          - R-GUI: drag and drop this file into an R-GUI window
##        The script will loop over all gct files in gct directory and run ssGSEA on each file 
##        separately.
################################################################################################################
setwd("");
rm(list=ls());
require("pacman");
require("cmapR");
script.dir = "";
## ##########################################################
##  define parameters below:
## ##########################################################
## ssGSEA / PTM-SEA parameters
sample.norm.type    = "none"              ## "rank", "log", "log.rank", "none" 
weight              = 0                ## value between 0 (no weighting) and 1 (actual data counts)
statistic           = "area.under.RES"    ## "Kolmogorov-Smirnov"
output.score.type   = "NES"               ## 'ES' or 'NES'
nperm               = 1e3                ## No. of permutations
min.overlap         = 3                  ## minimal overlap between gene set and data
correl.type         = "z.score"           ## 'rank', 'z.score', 'symm.rank'
par                 = T                   ## use 'doParallel' package?
spare.cores         = 3                   ## No. of cores to leave idle
export.signat.gct   = F                   ## if TRUE gene set GCT files will be exported 
extended.output     = F                   ## if TRUE the GCT files will contain stats on gene set overlaps etc.   
## #####################################################################
##   end paramaters
## - in a perfect world users don't have to worry about the stuff below...
## #####################################################################

## #################################
ORGANISEM = "";
GENESET = "genesetsYeast.gmt";
OUTDIR = "Yeast";
## directory with gct files
gct.dir <- paste0(script.dir, "/gct/", ORGANISEM);
## directory to write output
out.dir <- paste0(script.dir, "/out/", ORGANISEM);
## MSigDB
gene.set.databases <- paste0(script.dir, "/gmt/", GENESET);
## ######################################################################
##                          START
## ######################################################################
source(paste0(script.dir,'/src/ssGSEA2.0.R'));
## #############################################
## prepare output folder
setwd(out.dir)
date.str <- paste(OUTDIR, sep='_');
dir.create(date.str);
setwd(date.str);
## #############################################
## import signature database
signat.all <- unlist(lapply(gene.set.databases, readLines));
signat.all <- strsplit(signat.all, '\t');
names(signat.all) <- sapply(signat.all, function(x)x[1]);
signat.all <- lapply(signat.all, function(x) x[-c(1,2)]);

## save parameters used for ssGSEA
param.str = c(
    paste('##', Sys.time()),
    paste('gct.directory:', gct.dir, sep='\t'),
    paste('output.directory:', out.dir, sep='\t'),
    paste('gene.set.database:',gene.set.databases, sep='\t'),
    paste('sample.norm.type:', sample.norm.type, sep='\t'),
    paste('weight:', weight, sep='\t'),
    paste('statistic:', statistic, sep='\t'),
    paste('output.score.type', output.score.type, sep='\t'),
    paste('nperm:', nperm, sep='\t'),
    paste('min.overlap:', min.overlap, sep='\t'),
    paste('correl.type:', correl.type, sep='\t'),
    paste('run.parallel:', par, sep='\t')
   )
writeLines(param.str, con='parameters.txt')

## identify all gct files
gct.files <- dir(gct.dir, pattern='\\.gct$', full.names=T)
names(gct.files) <- paste(  sub('\\.gct$', '', sub('.*/','', gct.files)), 'ssGSEA', sep='_' )

## #####################################
## loop over gct files and run ssGSEA
for(i in names(gct.files)){
    ## create sub folders if more than one gct file was found
    if(length(gct.files) > 1){
        subdir=sub(',|\\.|:|;|/', '_', i)
        dir.create(subdir)
        setwd(subdir)
    }

    ## ########################################
    ## ssGSEA

    ## input data set
    input.ds <- gct.files[i]

    cat('Running ssSGEA on:', sub('.*/', '', input.ds), '\n\n')

    ## run ssGSEA
    gsea.res <- ssGSEA2(input.ds, gene.set.databases=gene.set.databases, sample.norm.type=sample.norm.type, weight=weight,statistic=statistic, output.score.type = output.score.type, nperm  = nperm, min.overlap  = min.overlap, correl.type = correl.type, output.prefix = paste(i), par=par, 
                        spare.cores=spare.cores, param.file=F, export.signat.gct = export.signat.gct, extended.output = extended.output )

    ## save object
    #save(gsea.res, file=paste(i, '.RData', sep=''))

    if(length(gct.files) > 1)
        setwd('..')

}






