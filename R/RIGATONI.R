
#' Create an upstream gene list, and a downstream gene list using STRING
#'
#' @param gene A character string naming the gene of interest
#' @param sdb The list of protein interactions (provided)
#'
#' @return A list of two vectors with the upstream and downstream results
#' @export
makeGeneList <- function(gene, sdb = sdb) {
  #sbd is the string database file
  #gene is the name of the gene of interest
  #filter sdb to only include rows which include the gene of interest
  sdb_goi <- sdb[sdb$item_id_a == gene |
                   sdb$item_id_b == gene, ]
  #remove rows where the acting gene (item_id_a) is absent
  sdb_goi = sdb_goi[!(is.na(sdb_goi$item_id_a)), ]
  #create the upstream database by searching for cases where the acting gene is gene a, the gene being acted on is the gene of interest, and the mode is expression.
  upstream = sdb_goi[sdb_goi$item_id_b == gene &
                       sdb_goi$a_is_acting == 't' &
                       sdb_goi$mode == 'expression', ]
  #create one version of downstream where the gene of interest is acting
  downstream_1 = sdb_goi[sdb_goi$item_id_a == gene &
                           sdb_goi$a_is_acting == 't', ]
  #create another downstream where the gene of interest is not necesarily acting, but expression is not being effected
  downstream_2 = sdb_goi[sdb_goi$item_id_b == gene &
                           sdb_goi$a_is_acting == 't' &
                           sdb_goi$mode != 'expression', ]
  #combine the two downstream gene lists
  downstream = c(downstream_2$item_id_a, downstream_1$item_id_b)
  #extract the gene list from the upstream data
  upstream = upstream$item_id_a
  #add both the upstream and downstream lists to a master list l
  l = list(upstream, downstream)
  #name both entries
  names(l) = c('upstream', 'downstream')
  #return the master list
  return(l)
}

#' Create two regression functions with the vectors from makeGeneList
#'
#' @param gene_list_ppi makeGeneList output
#' @param gene A character string naming the gene of interest
#' @param ControlRNA wild type control RNA seq counts where columns are samples and rows are genes
#'
#' @return A list of two linear regression models
#' @export
getRegression <- function(gene_list_ppi, gene, ControlRNA){
  #check that the length of the upstream gene list is greater than 0
  if (length(gene_list_ppi$upstream) > 0){
    #if it is greater than 0, build the regression model for the upstream list
    RegressionUpstream <- buildRegression(ControlRNA, gene, gene_list_ppi$upstream)
  } else {
    #if not, call the regression model "skip"
    RegressionUpstream = 'Skip'
  }
  #check that the length of the downstream gene list is greater than 0
  if (length(gene_list_ppi$downstream) > 0){
    #if it is greater than 0, build the regression model for the downstream list
    RegressionDownstream <- buildRegression(ControlRNA, gene, gene_list_ppi$downstream)
  } else {
    #if not, call the regression model "skip"
    RegressionDownstream = 'Skip'
  }
  #store each regression model in a master list
  l = list(RegressionUpstream, RegressionDownstream)
  #name each model
  names(l) = c('RegressionUpstream', 'RegressionDownstream')
  #return the list
  return(l)
}

#' Create one regression functions with one vector of predictor genes
#'
#' @param gene_list A character vector of genes to use as predictors
#' @param gene A character string naming the gene of interest
#' @param ControlRNA wild type control RNA seq counts where columns are samples and rows are genes
#'
#' @return A linear regression model
#' @export
buildRegression <- function(ControlRNA, gene, gene_list){
  #create a vector of the counts of the gene of interest called cd
  cd = t(ControlRNA[rownames(ControlRNA) == gene, ])
  #name the entries after the colnames of the control RNA
  rownames(cd) = colnames(ControlRNA)
  #turn the vector into a dataframe
  cd = as.data.frame(cd)
  #create the x.vars as a dataframe with genes as columns and sample names as rows
  #the genes in this dataframe should only be in the gene_list provided
  x.vars = as.data.frame(t(ControlRNA[rownames(ControlRNA) %in% gene_list, ]))
  #combine the x.vars dataframe with the gene of interest counts
  x.vars = cbind(x.vars, cd[, 1])
  #rename the gene of interest column GOI
  colnames(x.vars)[ncol(x.vars)] = 'GOI'
  #replace any NA entries with 0
  x.vars[is.na(x.vars)] = 0
  #create the model using the GOI column as a response and all other columns as predictors using a poisson distribution
  model = stats::glm(GOI ~ ., data = x.vars, family = stats::poisson())
  #return the model
  return(model)
}

#' Get prediction intervals of gene of interest counts for mutant samples
#'
#' @import ciTools
#' @param Regression getRegression output
#' @param gene A character string naming the gene of interest
#' @param MasterRNA mutated samples' RNA seq counts where columns are samples and rows are genes
#'
#' @return A dataframe with the true counts and upper and lower 95% prediction interval bounds
#' @export
mutantRegression <- function(Regression, MasterRNA, gene) {
  #generate the gene names in the upstream and downstream regressions
  if (Regression$RegressionUpstream != 'Skip'){
    names_up = names(Regression$RegressionUpstream$coefficients)
  } else {
    #if the upstream regression was skipped, make the names_up vector empty
    names_up = c()
  }
  if (Regression$RegressionDownstream != 'Skip'){
    names_down = names(Regression$RegressionDownstream$coefficients)
  } else {
    #if the downstream regression was skipped, make the names_up vector empty
    names_down = c()
  }
  #store the genes in upstream and downstream together in gene_list
  gene_list <- c(names_up, names_down)
  #store the regressions respectively
  upstream = Regression$RegressionUpstream
  downstream = Regression$RegressionDownstream
  #create a dataframe with the counts of the gene of interest (cd)
  cd = t(MasterRNA[rownames(MasterRNA) == gene, ])
  #save the sample names as rownames of the GOI dataframe (cd)
  rownames(cd) = colnames(MasterRNA)
  cd = as.data.frame(cd)
  #create the x.vars as a dataframe with genes as columns and sample names as rows
  #the genes in this dataframe should only be in the gene_list provided
  x.vars = as.data.frame(t(MasterRNA[rownames(MasterRNA) %in% gene_list, ]))
  #rename the gene of interest column GOI
  x.vars$GOI = cd[, 1]
  #create prediction intervals for each sample, for each regression model
  #alpha = .1 means you are creating a 95% prediction interval
  #nsims is the number of simulations of prediction to complete
  if (upstream != 'Skip'){
    fitted_up <- suppressWarnings(ciTools::add_pi(x.vars, upstream, alpha = 0.1, nsims = 20000))
  } else {
    #if no upstream model is present, skip it
    fitted_up <- 'Skip'
  }
  if (downstream != 'Skip'){
    fitted_down <- suppressWarnings(ciTools::add_pi(x.vars, downstream, alpha = 0.1, nsims = 20000))
  } else {
    #if no downstream model is present, skip it
    fitted_down <- 'Skip'
  }
  #put the fittered dataframes in to a list called l
  l = list(fitted_up, fitted_down)
  #for each dataframe in l
  l = lapply(l, function(f){
    if (f != 'Skip'){
      #as long as the entry of l is not "skip"
      #select only the columns listed below
      f = f[, colnames(f) %in% c('GOI', 'LPB0.05', 'UPB0.95')]
      #initialize a new vector called anno
      anno =c()
      for (x in 1:nrow(f)){
        #for each sample
        #if the true GOI count is outside the 95% prediction interval mark it F
        if (f$GOI[x] < f$LPB0.05[x] |
            f$GOI[x] > f$UPB0.95[x]){
          anno = c(anno, F)
        } else {
          anno = c(anno, T)
        }
      }
      #add a column to f called Annotation with the entries of anno
      f$Annotation = anno
      return(f)
    } else {
      #if the model was empty, skip it
      return('Skip')
    }
  })
  #rename the entries of l
  names(l) = c('upstream', 'downstream')
  #return l
  return(l)
}

#' Predict function of gene of interest within mutant samples
#'
#' @param Regression mutantRegression output
#' @param gene A character string naming the gene of interest
#' @param MasterRNA mutated samples' RNA seq counts where columns are samples and rows are genes
#'
#' @return A data.frame with sample names (taken from MasterRNA column names) and the predicted function of gene of interest
#' @export
getMutantFunction <- function(Regression, MasterRNA, gene){
  #first the function using mutantRegression function (above) to determine which samples are outside the 95% prediction interval
  Mutant_Regression <- suppressWarnings(mutantRegression(Regression, MasterRNA, gene))
  #remove the entries of the list that we should "skip"
  Mutant_Regression_1 = Mutant_Regression[Mutant_Regression != 'Skip']
  if (length(Mutant_Regression_1) > 1){
    #if both regression models are working, combine them into a master data frame called Anno
    Anno <- do.call(cbind, Mutant_Regression)
    Anno = as.data.frame(Anno)
  } else {
    #if at one regression model is skipped, renmae the remaining model "Anno"
    Anno = Mutant_Regression_1[[1]]
  }
  #initialize empty vector fun
  fun <- c()
  if (Mutant_Regression$downstream != 'Skip'){
    #if the downstream regression was able to be created continue below
    for (x in 1:nrow(Anno)){
      #for each sample
      #create new vector with the GOI counts
      downstream.GOI = ifelse((T %in% grepl('downstream', colnames(Anno))),
                              Anno$downstream.GOI[x],
                              Anno$GOI[x])
      #create new vector with lower prediction interval bound
      downstream.LPB0.05 = ifelse((T %in% grepl('downstream', colnames(Anno))),
                                  Anno$downstream.LPB0.05[x],
                                  Anno$LPB0.05[x])
      #create new vector with upper prediction interval bound
      downstream.UPB0.95 = ifelse((T %in% grepl('downstream', colnames(Anno))),
                                  Anno$downstream.UPB0.95[x],
                                  Anno$UPB0.95[x])
      if (downstream.GOI < downstream.LPB0.05){
        #if the true value is less than the lower bound, this is a LOF sample
        fun = c(fun, 'LOF')
      } else if (downstream.GOI > downstream.UPB0.95){
        #if the true value is greater than the upper bound, this is a GOF sample
        fun = c(fun, 'GOF')
      } else {
        #if the downstream predictions are within the boundaries, keep checking with the upstream samples
        fun = c(fun, 'keepchecking')
      }
    }
    #remove the solved samples, call them solved
    AnnoSolved <- Anno[which(fun != 'keepchecking'), ]
    #create new column with the function data
    AnnoSolved$Function = fun[which(fun != 'keepchecking')]
    #put the samples that need to be rechecked into their own dataframe
    AnnoKeep <- Anno[which(fun == 'keepchecking'), ]
  } else {
    #if there are no downstream values, annosovled should be empty
    AnnoSolved = NULL
    #annokeep should be the entirety of anno
    AnnoKeep = Anno
  }
  if (Mutant_Regression$upstream != 'Skip'){
    #if the upstream regression is present
    #initialize the empty fun vector
    fun = c()
    for (x in 1:nrow(AnnoKeep)){
      #for each sample
      #create new vector with the GOI counts
      upstream.GOI = ifelse((T %in% grepl('upstream', colnames(Anno))),
                            AnnoKeep$upstream.GOI[x],
                            AnnoKeep$GOI[x])
      #create new vector with lower prediction interval bound
      upstream.LPB0.05 = ifelse((T %in% grepl('upstream', colnames(Anno))),
                                AnnoKeep$upstream.LPB0.05[x],
                                AnnoKeep$LPB0.05[x])
      #create new vector with upper prediction interval bound
      upstream.UPB0.95 = ifelse((T %in% grepl('upstream', colnames(Anno))),
                                AnnoKeep$upstream.UPB0.95[x],
                                AnnoKeep$UPB0.95[x])
      if (upstream.GOI < upstream.LPB0.05){
        #If the true value is less than lowerbound, this is a LOF sample
        fun = c(fun, 'LOF')
      } else if (upstream.GOI > upstream.UPB0.95){
        #if the true value is greater than the upperbound, this is a GOF sample
        fun = c(fun, 'GOF')
      } else {
        #if we cannot decide, call it skipq
        fun = c(fun, "Skipq")
      }
    }
    #replace skipq with "unknown
    fun[fun == 'Skipq'] = 'Unknown'
    #create new column called Function with the vector
    AnnoKeep$Function = fun
  }
  if (!('Function' %in% colnames(AnnoKeep))){
    #if we had no upstream regression, initialize columnn Function with all entries "unknown"
    AnnoKeep$Function = 'Unknown'
  }
  #combine the solved and kept dataframes into Mutant_Function
  Mutant_Function <- rbind(AnnoSolved, AnnoKeep)
  #remove extraneous columns
  Mutant_Function <- as.data.frame(cbind(rownames(Mutant_Function), Mutant_Function$Function))
  #return mutant function
  return(Mutant_Function)
}

#' Predict immune phenotype within mutant samples
#'
#' @import ordinalForest
#' @import glmnet
#' @import DGEobj.utils
#' @param MasterRNA mutated samples' RNA seq counts where columns are samples and rows are genes
#' @param lasso_best ElasticNet output to select genes correllated with immune phenotype (provided)
#' @param m OrdinalForest model to predict immune phenotype (provided)
#' @param geneLen Dataframe of gene lengths (provided)
#'
#' @return A data.frame with sample names (taken from MasterRNA column names) and the predicted immune phenotype
#' @export
getImmuneProb <- function(MasterRNA, lasso_best = lasso_best, m = m , geneLen = geneLen){
  geneLen = geneLen[geneLen$Gene_Symbol %in% rownames(MasterRNA), ]
  MasterRNA = MasterRNA[rownames(MasterRNA) %in% geneLen$Gene_Symbol, ]
  geneLen = geneLen[match(rownames(MasterRNA), geneLen$Gene_Symbol), ]
  rnatmp = DGEobj.utils::convertCounts(as.matrix(MasterRNA), 'TPM', geneLen$Length)
  rnatmp = as.data.frame(rnatmp)
  genes = as.matrix(glmnet::coef.glmnet(lasso_best))
  genes = genes[!(genes[, 1] == 0), ]
  genes = names(genes)
  rna = rnatmp[rownames(rnatmp) %in% genes, ]
  col = rownames(rna)
  rna = apply(rna, 2, scale)
  rna = as.data.frame(t(rna))
  colnames(rna) = col
  pred <- stats::predict(m, newdata = rna, type = 'class')
  pred = cbind(rownames(rna), pred$classprobs[, 2])
  pred = as.data.frame(pred)
  colnames(pred) = c('Names', 'Prob')
  pred$Prob = as.numeric(as.character(pred$Prob))
  return(pred)
}

#' Predict immune phenotype within mutant samples
#'
#' @param gene A character string naming the gene of interest
#' @param ControlRNA wild type control RNA seq counts where columns are samples and rows are genes
#' @param MasterRNA mutated samples' RNA seq counts where columns are samples and rows are genes
#' @param lasso_best ElasticNet output to select genes correllated with immune phenotype (provided)
#' @param m OrdinalForest model to predict immune phenotype (provided)
#' @param geneLen Dataframe of gene lengths (provided)
#' @param sdb The list of protein interactions (provided)
#'
#' @return A data.frame with sample names (taken from MasterRNA column names), function prediction, and the predicted immune phenotype
#' @export
runRIGATONI <- function(gene, ControlRNA, MasterRNA,
                        lasso_best = lasso_best, m = m,
                        geneLen = geneLen,
                        sdb = sdb){
  if (any(class(ControlRNA) != 'data.frame' | class(MasterRNA) != 'data.frame')) stop('RNA files should be data.frame objects')
  if (rownames(ControlRNA) != rownames(MasterRNA)) stop('Gene names in ControlRNA are not the same as MasterRNA')
  if (F %in% (rownames(ControlRNA) %in% geneLen$Gene_Symbol)) stop('Genes in ControlRNA are not all present in geneLen file. \nPlease generate your own geneLen file')
  if (!(gene %in% rownames(ControlRNA))) stop('Gene of interest count data is not included in RNA provided')
  if (!(gene %in% geneLen$Gene_Symbol)) stop('Gene of interest is not included in the default gene lengths dataframe. \nPlease provide your own gene lengths to continue')
  message(paste0('Starting gene: ', gene))
  message('Making Initial Gene List')
  gene_list_ppi = makeGeneList(gene, sdb = sdb)
  if (length(unlist(gene_list_ppi)) > 0) stop('Sadly, we do not have enough protein connections in STRING to analyze you gene-of-interest using runRIGATONI. \nPlease investigate other RIGATONI functions to see if another analysis process would be more suitable.')
  message('Building regression model')
  Regression <- suppressWarnings(getRegression(gene_list_ppi, gene, ControlRNA))
  message('Predicting mutant function')
  Mutant_Function <- suppressWarnings(getMutantFunction(Regression, MasterRNA, gene))
  message('Predicting immune phenotype')
  ImmuneProb <- suppressMessages(getImmuneProb(MasterRNA, lasso_best = lasso_best, m = m, geneLen = geneLen))
  message('Getting final output')
  out = cbind(colnames(MasterRNA), Mutant_Function[, 2], ImmuneProb[, 2])
  colnames(out) = c('SampleID', 'Function', 'ImmunePhenotype')
  out = as.data.frame(out)
  message('Done!')
  return(out)
}
