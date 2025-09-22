#' Create an upstream gene list, and a downstream gene list using STRING
#'
#' @param gene A character string naming the gene of interest
#' @param string The list of protein interactions (provided)
#'
#' @return A list of two vectors with the upstream and downstream results
#' @export
makeGeneList <- function(gene, string = sdb) {
  #sbd is the string database file
  #gene is the name of the gene of interest
  if (!(gene %in% c(string$item_id_a, string$item_id_b))) stop('Sadly, we do not have enough protein connections in STRING to analyze you gene-of-interest using runRIGATONI. \nPlease build your own gene list using prior knowledge or literature review.')
  #filter sdb to only include rows which include the gene of interest
  sdb_goi <- string[string$item_id_a == gene |
                      string$item_id_b == gene, ]
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
  downstream = downstream[downstream != '']
  #extract the gene list from the upstream data
  upstream = upstream$item_id_a
  upstream = upstream[upstream != '']
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
  if (length(Regression$RegressionUpstream) > 1 ) {
    names_up = names(Regression$RegressionUpstream$coefficients)
  } else {
    #if the upstream regression was skipped, make the names_up vector empty
    names_up = c()
  }
  if (length(Regression$RegressionDownstream) > 1 ) {
    names_down = names(Regression$RegressionDownstream$coefficients)
  } else {
    #if the downstream regression was skipped, make the names_up vector empty
    names_down = c()
  }
  #store the genes in upstream and downstream together in gene_list
  gene_list <- c(names_up, names_down)
  rownames(MasterRNA) = gsub("\\-", "\\.", rownames(MasterRNA))
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
  if (length(upstream) > 1) {
    fitted_up <- suppressWarnings(ciTools::add_pi(x.vars, upstream, alpha = 0.1, nsims = 20000))
  } else {
    #if no upstream model is present, skip it
    fitted_up <- 'Skip'
  }
  if (length(downstream) > 1) {
    fitted_down <- suppressWarnings(ciTools::add_pi(x.vars, downstream, alpha = 0.1, nsims = 20000))
  } else {
    #if no downstream model is present, skip it
    fitted_down <- 'Skip'
  }
  #put the fittered dataframes in to a list called l
  l = list(fitted_up, fitted_down)
  #for each dataframe in l
  l = lapply(l, function(f){
    if (length(f) > 1){
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
  if (length(Mutant_Regression$downstream) > 1){
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
  if (length(Mutant_Regression$upstream) > 1){
    #if the upstream regression is present
    #initialize the empty fun vector
    fun = c()
    #as long as the downstream didn't solve everything
    if (nrow(AnnoKeep > 0)) {
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
    }
    #replace skipq with "unknown
    fun[fun == 'Skipq'] = 'Unknown'
    #create new column called Function with the vector
    AnnoKeep$Function = fun
  }
  if (!('Function' %in% colnames(AnnoKeep))){
    #if we had no upstream regression, initialize columnn Function with all entries "unknown"
    AnnoKeep$Function = rep('Unknown', nrow(AnnoKeep))
  }
  #combine the solved and kept dataframes into Mutant_Function
  Mutant_Function <- rbind(AnnoSolved, AnnoKeep)
  #remove extraneous columns
  Mutant_Function <- as.data.frame(cbind(rownames(Mutant_Function), Mutant_Function$Function))
  #return mutant function
  return(Mutant_Function)
}

#' Predict immune phenotype
#'
#' @param MasterRNA samples' RNA seq transcripts per million where columns are samples and rows are genes
#' @param model_l Model to predict Low immune advantage immune phenotype (provided)
#' @param model_h Model to predict High immune advantage immune phenotype (provided)
#' @param use_genes genes to predict Low immune advantage immune phenotype (provided)
#' @importFrom stats predict
#' @importFrom pamr pamr.predict
#'
#' @return A data.frame with sample names (taken from MasterRNA column names) and the predicted immune phenotype
#' @export
getImmuneProb <- function(MasterRNA, model_l = mod_l, model_h = mod_h, use_genes = genes) {
  
  #scale the RNA and transpose it
  MasterRNA = t(scale(MasterRNA))
  
  # get class probabilities
  predictions = data.frame(High = pamr.predict(model_h, t(MasterRNA[, use_genes$High]), threshold = 0, type = 'posterior')[, 1],
                           Low = predict(model_l, MasterRNA[, use_genes$Low], type = 'prob')[, 1])
  
  # pick the right prediction for each row
  pred = unlist(lapply(1:nrow(predictions), function(x){
    
    # the low model was better performing than the high model, so pick low first
    if (predictions[x, 2] > .5) {
      
      return('Low')
      
    } else if (predictions[x, 1] > .5) {
      
      return('High')
      
    } else {
      
      # if a sample isn't High or Low, call it "Ind" or indeterminant
      
      return('Ind')
      
    }
  }))
  
  #Add sample names to the prediction
  pred = cbind(rownames(MasterRNA), pred)
  
  #create data frame
  pred = as.data.frame(pred)
  
  #add column names
  colnames(pred) = c('Names', 'Prob')
  
  #return data
  return(pred)
  
}

#' Evaluate variant harbored by mutant samples
#'
#' @param Function A data.frame returned by getMutantFunction
#' @param ImmuneProb A data.frame returned by getImmuneProb
#' @param ControlRNA wild type control RNA seq counts where columns are samples and rows are genes
#' @importFrom stats prop.test
#' 
#' @return A list with the elements Function (function annotation of variant), p.val.func (p value of function annotation), ImmunePhenotype (immune phenotype annotation of variant), and p.val.immune (p value of immune phenotype annotation)
#' @export
evaluateMutants <- function(Function, ImmuneProb, ControlRNA){
  con = suppressMessages(getImmuneProb(ControlRNA, model_l = mod_l, model_h = mod_h, use_genes = genes))
  nGOF = nrow(Function[Function[, 2] == 'GOF', ])
  nLOF = nrow(Function[Function[, 2] == 'LOF', ])
  nHigh = nrow(ImmuneProb[ImmuneProb[, 2] == 'High', ])
  nLow = nrow(ImmuneProb[ImmuneProb[, 2] == 'Low', ])
  conPropH = nrow(con[con[, 2] == 'High', ])/nrow(con)
  conPropL = nrow(con[con[, 2] == 'Low', ])/nrow(con)
  pf = prop.test(nGOF, sum(nGOF, nLOF), .5)
  Fun = ifelse(pf$p.value < .05, ifelse(nGOF > nLOF, 'GOF', 'LOF'), 'Unknown')
  pih = prop.test(nHigh, nrow(ImmuneProb), conPropH, 'greater')
  pil = prop.test(nLow, nrow(ImmuneProb), conPropL, 'greater')
  Imm = ifelse(pih$p.value < .05, 'High', ifelse(pil$p.value < .05, 'Low', 'Unknown'))
  pi = ifelse(Imm == 'High', pih$p.value, ifelse(Imm == 'Low', pil$p.value, min(pil$p.value, pih$p.value)))
  out = list(Fun, pf, Imm, pi)
  names(out) = c('Function', 'p.val.func', 'ImmunePhenotype', 'p.val.immune')
  return(out)
}

#' Predict immune phenotype within mutant samples
#'
#' @param gene A character string naming the gene of interest
#' @param ControlRNA wild type control RNA seq counts where columns are samples and rows are genes
#' @param MasterRNA mutated samples' RNA seq counts where columns are samples and rows are genes
#' @param model_l Model to predict Low immune advantage immune phenotype (provided)
#' @param model_h Model to predict High immune advantage immune phenotype (provided)
#' @param use_genes genes to predict Low immune advantage immune phenotype (provided)
#' @param string The list of protein interactions (default provided)
#' @param geneList list of upstream and downstream genes to use to predict mutant function. Default value is NULL and geneList will be built using STRING database.
#'
#' @return A data.frame with sample names (taken from MasterRNA column names), function prediction, and the predicted immune phenotype
#' @export
runRIGATonI <- function(gene, ControlRNA, MasterRNA,
                        model_l = mod_l, 
                        model_h = mod_h, 
                        use_genes = genes,
                        string = sdb,
                        geneList = NULL){
  if (any(class(ControlRNA) != 'data.frame', class(MasterRNA) != 'data.frame')) stop('RNA files should be data.frame objects')
  if (any(!(rownames(ControlRNA) == rownames(MasterRNA)))) stop('Gene names in ControlRNA are not the same as MasterRNA')
  if (!(gene %in% rownames(ControlRNA))) stop('Gene of interest TPM data is not included in RNA provided')
  if (any(!(unlist(genes) %in% rownames(MasterRNA)))) stop('The genes required for immune phenotype calculation are not all present in the data provided. \nPlease check that all your gene names are uppercase and the rownames of your RNA TPM data.')
  message(paste0('Starting gene: ', gene))
  if (is.null(geneList) == T ) {
    message('Making Initial Gene List')
    gene_list_ppi = makeGeneList(gene, string = sdb)
    message('Filtering Gene List')
    gene_list_ppi = lapply(gene_list_ppi, function(x){
      out = x[x %in% rownames(ControlRNA)]
      return(out)
    })
    if (length(unlist(gene_list_ppi)) < 0) stop('Sadly, we do not have enough protein connections in STRING to analyze you gene-of-interest using runRIGATONI. \nPlease build your own gene list using prior knowledge or literature review.')
  } else {
    gene_list_ppi = geneList
    message('Filtering Gene List')
    gene_list_ppi = lapply(gene_list_ppi, function(x){
      out = x[x %in% rownames(ControlRNA)]
      return(out)
    })
    if (length(unlist(gene_list_ppi)) > 0) stop('Your gene list is either empty or does not contain genes in your RNA TPM data. \nPlease re-run RIGATONI with default settings and allow us to make the geneList using STRING.')
  }
  message('Building regression model')
  Regression <- suppressWarnings(getRegression(gene_list_ppi, gene, ControlRNA))
  message('Predicting mutant function')
  Mutant_Function <- suppressWarnings(getMutantFunction(Regression, MasterRNA, gene))
  message('Predicting immune phenotype')
  ImmuneProb <- suppressMessages(getImmuneProb(MasterRNA, model_l = mod_l, model_h = mod_h, use_genes = genes))
  message('Getting final output')
  out = cbind(colnames(MasterRNA), Mutant_Function[, 2], ImmuneProb[, 2])
  colnames(out) = c('SampleID', 'Function', 'ImmunePhenotype')
  out = as.data.frame(out)
  message('Calculating Alteration Function and Phenotype')
  muts = evaluateMutants(Mutant_Function, ImmuneProb, ControlRNA)
  message(cat(paste0(
    'Function of Variant: ', muts$Function, ", p = ", muts$p.val.func$p.value, "\n",
    "Immune Phenotype of Variant: ", muts$ImmunePhenotype, ", p = ", muts$p.val.immune
  )))
  message('Done!')
  return(out)
}
