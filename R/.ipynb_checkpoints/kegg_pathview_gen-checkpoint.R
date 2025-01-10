## Generate KEGG pathviews from SCODA result
## Seokhyun Yoon (syoon@dku.edu), MLBI@DKU, Jan 03, 2025

# suppressPackageStartupMessages(library(stringr))
# suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(filesstrings))
# suppressPackageStartupMessages(library(pathview))
# suppressPackageStartupMessages(library(gage))
# suppressPackageStartupMessages(library(gageData))
# suppressPackageStartupMessages(library(org.Hs.eg.db))
# suppressPackageStartupMessages(library(org.Mm.eg.db))
# suppressPackageStartupMessages(library(reticulate))
# suppressPackageStartupMessages(library(anndata))

add_entrez <- function( df.deg.res, species )
{
    if(species == 'hsa'){
        org.db <- org.Hs.eg.db
    } else{
        org.db <- org.Mm.eg.db
    }

    my.symbols <- rownames(df.deg.res)
    suppressMessages( Convert_result <- biomaRt::select(org.db, keys = my.symbols,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL") )
    Convert_result <- distinct(Convert_result, SYMBOL, .keep_all= TRUE)

    df.deg.res <- cbind(df.deg.res, Convert_result)
    return(df.deg.res)
}

get_species <- function( species )
{
    if( tolower(species) %in% c('hs', 'human') ){ species = 'hsa' }
    else if( tolower(species) %in% c('mm', 'mouse')){ species = 'mmu' }
    return(species)
}

get_fold_changes <- function( df.deg.res, species, pval.cutoff = 0.01 )
{
    if( tolower(species) %in% c('hs', 'human') ){ species = 'hsa' }
    else if( tolower(species) %in% c('mm', 'mouse')){ species = 'mmu' }

    if( !("ENTREZID" %in% colnames(df.deg.res)) )
    {
        df.deg.res <- add_entrez( df.deg.res, species )
    }

    b <- (df.deg.res$pval <= pval.cutoff)
    foldchanges = df.deg.res$log2_FC[b]
    names(foldchanges) = df.deg.res$ENTREZID[b]

    return(foldchanges)
}


#' get_all_fold_changes
#'
#' This function retreives fold changes from the SCODA DEG results.
#'
#' @param lst.deg.res DEG results from SCODA (a named list of data frames)
#' @param species Either 'human' or 'mouse'
#' @param pval.cutoff p-value cutoff to filter DEGs
#' @return A named list of named list of fold changes. Names of fold changes are ENTREZ gene code
#' @export
get_all_fold_changes <- function( lst.deg.res, species, pval.cutoff = 0.01 )
{
    cat(sprintf('Getting fold changes .. \n'))
    flush.console()

    lst.fc.res <- c()
    for( j in 1:length(lst.deg.res) )
    {
        cat(sprintf('%30s: ', names(lst.deg.res)[j]))
        flush.console()

        lst.deg.t <- lst.deg.res[[j]]
        lst.fc.t <- list()
        for( i in 1:length(lst.deg.t) )
        {
            foldchanges <- get_fold_changes(lst.deg.t[[i]], species, pval.cutoff )
            lst.fc.t[[i]] <- foldchanges
            if( i == length(lst.deg.t) )
            {
                cat(sprintf('%s(%d)\n', names(lst.deg.t)[i], sum(foldchanges != 0)))
            } else{
                cat(sprintf('%s(%d), ', names(lst.deg.t)[i], sum(foldchanges != 0)))
            }
        }
        names(lst.fc.t) <- names(lst.deg.t)
        lst.fc.res[[j]] <- lst.fc.t 
    }
    names(lst.fc.res) <- names(lst.deg.res)
    cat(sprintf('Getting fold changes .. done. \n'))
    return( lst.fc.res )
}


get_kegg_pathways_info <- function( species )
{
    kegg.db <- kegg.gsets(species = species, id.type = "kegg", check.new = TRUE)
    kegg_pw_id_name <- names(kegg.db$kg.sets)

    kegg_pw_id <- vapply(strsplit(kegg_pw_id_name, " ", fixed = TRUE), "[", "", 1)
    kegg_pw_name <- as.vector(vapply(kegg_pw_id_name, function(x) substring(x, 10), character(1)))

    df_pw <- data.frame( pw_id = kegg_pw_id, pw_name = kegg_pw_name, pw_id_name = kegg_pw_id_name )
    rownames(df_pw) <- (kegg_pw_id_name)

    pathways_kegg_entrez <- list()
    for( i in 1:length(kegg_pw_id_name) )
    {
        pathways_kegg_entrez[[i]] <- kegg.db$kg.sets[[kegg_pw_id_name[i]]]
    }
    names(pathways_kegg_entrez) <- kegg_pw_id_name

    return( list(id_name_map = df_pw, gene.sets = pathways_kegg_entrez ) )
}


convert_gene_symbols_in_Pathways_DB_to_entrez <- function(pathways_used, species)
{
    if( tolower(species) %in% c('hs', 'human') ){ species = 'hsa' }
    else if( tolower(species) %in% c('mm', 'mouse')){ species = 'mmu' }

    if(species == 'hsa'){
        org.db <- org.Hs.eg.db
    } else{
        org.db <- org.Mm.eg.db
    }

    pathways_names_used <- names(pathways_used)
    ## For each pathway, convert hugo symbols into ENTREZ symbol
    pathways_used_entrez <- list()
    cat(sprintf('Converting Pathways DB .. \r'))
    for( i in 1:length(pathways_names_used) )
    {
        cat(sprintf('Converting Pathways DB .. %3d/%3d   \r', i, length(pathways_names_used)) ) #, pathways_names_used[i]))
        flush.console()

        hugo.symbols <- pathways_used[[pathways_names_used[i]]]
        if( length(hugo.symbols) == 1 ){ hugo.symbols <- hugo.symbols[[1]] }

        {
            s <- hugo.symbols[[1]][length(hugo.symbols[[1]])]
            if(  substr( s, str_length(s), str_length(s) ) == '\n' )
            {
                hugo.symbols[[1]][length(hugo.symbols[[1]])] <- substr( s, 1, str_length(s)-1 )
            }
        }

        suppressMessages( Convert_result <- biomaRt::select(org.db, keys = hugo.symbols,
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL") )
        Convert_result <- distinct(Convert_result, SYMBOL, .keep_all= TRUE)
        pathways_used_entrez[[i]] <- Convert_result[,"ENTREZID"]
        wh <- which(is.na( pathways_used_entrez[[i]] ))
        pathways_used_entrez[[i]] <- pathways_used_entrez[[i]][-wh]
    }
    names(pathways_used_entrez) <- pathways_names_used
    cat(sprintf('Converting Pathways DB .. done.        \n'))

    return(pathways_used_entrez)
}


## For each pathway, find the best-match KEGG pathways
build_pathway_map <- function(pathways_used_entrez, pathways_kegg_entrez, o_th = 0.85)
{
    # o_th <- 0.85
    pathways_map <- c()
    names_map <- c()

    cnt <- 0
    for( i in 1:length(pathways_used_entrez) )
    {
        ov <- c()
        or <- c()
        for( j in 1:length(pathways_kegg_entrez) )
        {
            pwi <- intersect( pathways_used_entrez[[i]], pathways_kegg_entrez[[j]] )
            ov[j] <- length(pwi)
            or[j] <- length(pwi) /min(length(pathways_used_entrez[[i]]), length(pathways_kegg_entrez[[j]]))
        }
        w <- which.max(ov)
        if(or[w] >= o_th)
        {
            cnt <- cnt + 1
            pathways_map[cnt] <- names(pathways_kegg_entrez)[w]
            names_map[cnt] <- names(pathways_used_entrez)[i]
        }
    }
    names(pathways_map) <- names_map

    return(pathways_map)
}


select_valid_pathways_from_gsa_result <- function( df.gsa, pathways_map, df_kegg_pw_map,
                                                   p.val.cutoff = 0.01, verbose = FALSE )
{
    p.val.cutoff.scoda <- p.val.cutoff
    pv.col <- 'pval'

    n_before <- dim(df.gsa)[1]
    b <- df.gsa[,pv.col] <= p.val.cutoff.scoda
    df.gsa <- df.gsa[b,]
    pw_detected_scoda <- ( pathways_map[rownames(df.gsa)] )
    n_after <- dim(df.gsa)[1]
    pw_detected_scoda <- unique(pw_detected_scoda)

    if( verbose )
    {
        cat(sprintf('# GSA items: %d -> %d -> %d -> ', n_before, n_after, length(pw_detected_scoda) ))
    }

    kegg_pw_id_name <- rownames(df_kegg_pw_map)
    pw_scoda <- sort(pw_detected_scoda)
    pw_kegg <- sort((kegg_pw_id_name))
    pw_common <- intersect(pw_scoda, pw_kegg)
    pw_scoda_only <- setdiff(pw_scoda, pw_common)
    pw_kegg_only <- setdiff(pw_kegg, pw_common)

    n_common <- length(pw_common)
    n_scoda <- length(pw_scoda)
    n_scoda_only <- length(pw_scoda_only)
    n_kegg <- length(pw_kegg)

    if( verbose )
    {
        cat(sprintf('%d intersection %d -> %d \n', n_scoda, n_kegg, n_common ))
    }
    df_pw_sel <- df_kegg_pw_map[pw_common, ]

    return(df_pw_sel)
}


#' get_pathways_map
#'
#' This function retreives mapping between the KEGG pathways and the pathways DB used. 
#' Using the gene sets defined in the two pathways DBs, it selects, for each pathway in the SCODA Gene Set Analysis result, the KEGG pathway with the best match. 
#'
#' @param pathways_used a named list of gene sets used for the SCODA Gene Set (Enrichment) analysis. It is stored in the SCODA result file
#' @param species 'human' or 'mouse'
#' @param min_overlap a number between 0 and 1. Used to filter the matched KEGG pathways. The overlap is defined by the number of genes common in the two pathways divided by the smallest gene set size
#' @return a data frame containing the mapping between the KEGG pathway and the pathway used in the SCODA
#' @export
get_pathways_map <- function(pathways_used, species, min_overlap = 0.85 )
{
    if( tolower(species) %in% c('hs', 'human') ){ species = 'hsa' }
    else if( tolower(species) %in% c('mm', 'mouse')){ species = 'mmu' }

    pathways_used_entrez <- convert_gene_symbols_in_Pathways_DB_to_entrez( pathways_used, species )

    lst <- get_kegg_pathways_info( species = species )
    df_kegg_pw_map <- lst$id_name_map
    pathways_kegg_entrez <- lst$gene.sets

    pathways_map <- build_pathway_map( pathways_used_entrez, pathways_kegg_entrez, o_th = min_overlap )

    pws <- names(pathways_map)

    pw_sel <- c()
    cnt <- 0
    for( pw in pws )
    {
        if( pathways_map[[pw]] %in% rownames(df_kegg_pw_map) )
        {
            cnt <- cnt + 1
            pw_sel[cnt] <- pw
        }
    }

    df_kegg_pw_map_sel <- df_kegg_pw_map[ pathways_map[pw_sel], ]
    df_kegg_pw_map_sel[,'pw_name_used'] <- pw_sel

    return(df_kegg_pw_map_sel)
    # return( list(kegg_pw_map = df_kegg_pw_map, pathways_map = pathways_map) )
}


#' save_kegg_pathviews
#'
#' This function saves Kthe EGG pathview images colored by log-fold-changes of the pathways identified in the SCODA gene set (enrichment) analysis result for a specific cell-type. 
#'
#' @param target_cell target cell-type to retreives its DEGs on the KEGG pathview images. Should be one of the name list of lst.gsa.all.
#' @param lst.gsa.all Gene Set (Enrichment) Analysis result retreived from the SCODA result file. 
#' @param lst.fcs.all A named list of log-fold-changes obtained by the function 'get_all_fold_changes'. 
#' @param df_pathways_map A data frame containing the mapping between the KEGG pathway and the pathway used in the SCODA. Use the output of the function 'get_pathways_map'. 
#' @param species 'human' or 'mouse'. 
#' @param gsa.p.val.cutoff p-value to filter the Gene Set (Enrichment) analysis results in the SCODA result file. 
#' @return The name of the folder where KEGG pathview images are stored. It is 'KEGG_pathview_(cell type)', where (cell type) is the 'target_cell'.
#' @export
save_kegg_pathviews <- function( target_cell, lst.gsa.all, lst.fcs.all, df_pathways_map, 
                                 species, gsa.p.val.cutoff = 0.01 )
{
    if( tolower(species) %in% c('hs', 'human') ){ species = 'hsa' }
    else if( tolower(species) %in% c('mm', 'mouse')){ species = 'mmu' }

    lst.df.gsa <- lst.gsa.all[[target_cell]]
    fcs.entrez <- lst.fcs.all[[target_cell]]

    pathways_map <- rownames(df_pathways_map)
    names(pathways_map) <- df_pathways_map[,'pw_name_used']

    df_kegg_pw_map <- df_pathways_map

    dir_to_save <- sprintf('KEGG_pathview_%s', target_cell )

    if( !file.exists(dir_to_save) ){
        dir.create(dir_to_save)
    } else{
        f.lst <- list.files(dir_to_save)
        for(i in 1:length(f.lst))
        {
            f <- paste0(dir_to_save, '/', f.lst[i])
            if( file.exists(f) ){ file.remove(f) }
        }
    }
    items <- names(lst.df.gsa)

    for( j in 1:length(items) )
    {
        item <- items[j]
        foldchanges <- fcs.entrez[[item]]

        df.gsa <- lst.df.gsa[[item]]
        if( dim(df.gsa)[1] > 0 )
        {
            df_pw_sel <- select_valid_pathways_from_gsa_result( df.gsa, pathways_map, df_kegg_pw_map,
                                                                p.val.cutoff = gsa.p.val.cutoff )

            pw_id_sel <- df_pw_sel$pw_id
            pw_name_sel <- df_pw_sel$pw_name

            s_max = 0
            for( i in 1:length(pw_name_sel)){ s_max <- max( s_max, str_length(pw_name_sel[i]) ) }

            for( i in 1:length(pw_name_sel))
            {
                pid <- pw_id_sel[i]
                pname <- pw_name_sel[i]
                pname <- str_replace( pname, '/', '_' )

                s_suffix = ''
                if( str_length(pname) < s_max )
                {
                    for( k in 1:(s_max - str_length(pname)) )
                    {
                        s_suffix <- paste0(s_suffix, ' ')
                    }
                }

                cat(sprintf('%30s: %d/%d - %d/%d - %s%s \r', target_cell, j, length(items),
                              i, length(pw_name_sel), pname, s_suffix ))
                flush.console()

                suppressMessages( pv.out <- pathview(gene.data = foldchanges, pathway.id=pid,
                          species=species, kegg.dir = dir_to_save, # out.suffix = 'pos',
                          low = list(gene = "turquoise", cpd = "blue"),
                          mid = list(gene = "gray", cpd = "gray"),
                          high = list(gene = "gold", cpd = "yellow"),
                          kegg.native = TRUE, same.layer = FALSE, min.nnodes = 5))

                if( is.list(pv.out) )
                {
                    file_out <- paste0(pname, '_', item, '.png')
                    file.rename(paste0(pid,'.pathview.png'), file_out)
                    suppressMessages( file.move(file_out, dir_to_save) )
                    if( file.exists(paste0(dir_to_save, '/', pid,'.png')) )
                    {
                        file.remove(paste0(dir_to_save, '/', pid,'.png'))
                    }
                    if( file.exists(paste0(dir_to_save, '/', pid,'.xml')) )
                    {
                        file.remove(paste0(dir_to_save, '/', pid,'.xml'))
                    }
                }
            }
        }
    }
    cat(sprintf('%30s: %d/%d - %d/%d - %s%s \n', target_cell, j, length(items),
                  i, length(pw_name_sel), pname, s_suffix ))
    return(dir_to_save)
}
