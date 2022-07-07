DMRforPairs <- function(classes_gene, classes_island, targetID, chr, 
    position, m.v, beta.v, min_n = 4, min_distance = 200, min_dM = 1.4, 
    recode = 1, sep = ";", method = "fdr", debug.v = FALSE, gs, do.parallel = 0)
{
    # merge gene/island classes
    classes <- merge_classes(classes_gene, classes_island, recode, sep)
    rownames(classes$pclass) <- rownames(m.v)
    rownames(classes$pclass_recoded) <- rownames(m.v)
    
    # find potential regions of interest
    regions <- regionfinder(targetID, chr, position, classes$pclass_recoded, 
        classes$no.pclass, classes$u_pclass, d = min_distance, m.v, beta.v, n_min = min_n, debug = debug.v, gs)
    
    # identify and test regions with differential methylation.
    message("Calculating statistics and testing differences between samples per region. Please be patient.")
    
    probes <- regions$perprobe
    m <- regions$valid.m
    b <- regions$valid.beta
    n <- dim(m)[2]
    dMth <- min_dM
    
    n <- dim(regions$valid.m)[2]
    ID <- data.frame(regions$boundaries$regionID)
    if (do.parallel == 0) {
        tested = apply(ID, 1, function(ID) testregion(ID, probes = probes, 
            m = m, b = b, n = n, dMth = dMth, do.format = FALSE))  #calculate statistics for each region
    }
    if (do.parallel == -1) {
        do.parallel = detectCores()
    }
    if (do.parallel > 0) {
        cl <- makeCluster(getOption("cl.cores", do.parallel))
        clusterExport(cl, c("probes", "m", "b", "dMth", "testregion", 
            "calc_stats", "n"), envir = environment())
        tested <- parApply(cl, ID, 1, function(ID) testregion(ID, 
            probes = probes, m = m, b = b, n = n, dMth = dMth, do.format = FALSE))  #calculate statistics for each region
        stopCluster(cl)
    }
    
    tested <- t(tested)
    i = combn(n, 2)  #testregion is called using apply without requesting proper formatting to minimize compute time. Formatting is done once on the complete output in the following lines.
    cn1 = paste("beta.median", colnames(b), sep = ".")
    cn2 = paste("m.median", colnames(m), sep = ".")
    cn3 = paste("median.delta.beta", colnames(m)[i[1, ]], "minus", colnames(m)[i[2, ]], sep = ".")
    cn4 = paste("median.delta.m", colnames(m)[i[1, ]], "minus", colnames(m)[i[2, ]], sep = ".")
    cn5 = paste("pairwise.p", colnames(m)[i[1, ]], "vs", colnames(m)[i[2, ]], sep = ".")
    colnames(tested) = c(cn1, cn2, cn3, cn4, cn5, "max.abs.median.delta", "p.value")
    tested <- as.data.frame(tested)
    tested$p.value.adjusted <- p.adjust(p = tested$p.value, method = method, 
        n = length(which(!is.na(tested$p.value))))  #calc adjusted p-values based on the actual nr of tests performed.
    tested <- cbind(regions$boundaries, tested)  #combine methylation level statistics with positional info of the regions
    
    return(list(classes = classes, regions = regions, tested = tested))
}

merge_classes <- function(refgene_class, island_class, recode = 1, sep = ";")
{
    message("Recoding annotation classes...")
    probes.pclass <- paste(refgene_class, island_class, sep = ";")
    probes.pclass.merged = probes.pclass
    if (length(recode) == 1) {
        if (recode == 0) {
            # cluster gene, tss and cpg island related probes
            a = c("Body", "5'UTR", "3'UTR", "1stExon", "TSS1500", 
                "TSS200", "Island", "N_Shelf", "N_Shore", "S_Shelf", 
                "S_Shore")
            recode = data.frame(a, a)
        } else if (recode == 1) {
            # retain detailed annotation scheme of illumina
            a = c("gene", "tss", "island")
            b = c("Body;5'UTR;3'UTR;1stExon", "TSS1500;TSS200", "Island;N_Shelf;N_Shore;S_Shelf;S_Shore")
            recode = data.frame(a, b)
        } else if (recode == 2) {
            # ignore classes
            probes.pclass.merged = rep("all", length(probes.pclass.merged)[1])
            a = c("all")
            b = c("all")
            recode = data.frame(a, b)
        }
    }
    
    valid.idx = c()
    counter = 0
    rc = data.frame(v1 = rep(NA, length(probes.pclass.merged)[1]), 
        stringsAsFactors = FALSE)
    
    for (i in 1:length(recode[, 1])) {
        counter = counter + 1
        ia = data.frame(strsplit(as.character(recode[i, 2]), ";"))
        for (j in 1:dim(ia)[1]) {
            message(paste(as.character(ia[j, 1]), "-->", as.character(recode[i, 
                1])))
            cur_rows = matchsymbol(probes.pclass.merged, as.character(ia[j, 
                1]), sep = ";")
            rc[cur_rows, counter] = as.character(recode[i, 1])
            valid.idx = union(valid.idx, cur_rows)
        }
    }
    
    rc_con <- do.call(paste, c(rc, sep = ";"))
    rc_con = as.matrix(rc_con)
    probes.no.pclass = seq(1, length(probes.pclass.merged), 1)[-unique(valid.idx)]
    
    u_pclass = unique(recode[, 1])
    
    return(list(pclass = as.matrix(probes.pclass.merged), pclass_recoded = rc_con, 
        no.pclass = probes.no.pclass, u_pclass = u_pclass))
}

regionfinder <- function(targetID, chr, position, pclass, r_excl, 
    u_pclass, d = 200, m.v, beta.v, n_min = 4, debug = FALSE, gs) # regionfinder
{
    gs = data.frame(gs)
    # exclude probes with no classification from m/beta/call/p-value
    # matrices
    if (length(r_excl) > 0) {
        probes = data.frame(1:length(targetID[-1 * r_excl]), targetID[-1 * 
            r_excl], chr[-1 * r_excl], position[-1 * r_excl], pclass[-1 * 
            r_excl])  #rowID = row index of probe in tables with unclassified probes excluded
        m.v = m.v[-1 * r_excl, ]
        beta.v = beta.v[-1 * r_excl, ]
        gs = gs[-1 * r_excl, ]
    } else {
        probes = data.frame(1:length(targetID), targetID, chr, position, 
            pclass)  #rowID = row index of probe in tables with unclassified probes excluded
    }
    colnames(probes) <- c("rowID", "probeID", "chr", "position", 
        "pClass")
    
    # prep the loops
    u_chr = unique(chr)
    u_chr = u_chr[order(u_chr)]
    n_chr = length(u_chr)
    n_pclass = length(u_pclass)
    regions = data.frame()
    regions2probes = matrix(NA, nrow = dim(probes[1]), ncol = n_pclass)
    regionID = 0
    
    if (debug == TRUE) {
            n_chr = 1
    }  #just for debug purposes, only run chr1
    
    # for each chromosome... (since the regionfinding process is
    # recursive, it was not passible to do this elegantly using # apply)
    for (curchr in 1:n_chr) {
        ptm <- proc.time()
        message(paste("Regionfinder: processing chr", u_chr[curchr], " (", curchr, "/", n_chr, ")", sep = ""))
        
        cur_probes = probes[which(as.character(probes$chr) == as.character(u_chr[curchr])), ]  #select probes on current chromosome
        cur_probes = cur_probes[order(cur_probes[4]), ]  #order on position (asc)
        
        for (cur_pclass in 1:n_pclass) {
            # for each class (per chromosome)... find regions based on
            # adjecent probes with same categories
            cur_pclass.idx = which(regexpr(u_pclass[cur_pclass], 
                cur_probes$pClass) > 0)  #find probes which have at least the current class in their annotation (content=index in cur_probes)
            cur_pclass.difference = diff(cur_pclass.idx)  #difference between idx of adjecent probes within the current class (>1 = probe(s) with other class in between)
            
            cur_pclass.boundaries.start_idx = cur_pclass.idx[which(cur_pclass.difference != 
                1) + 1]  #diff(1) = idx(2)-idx(1) --> diff(i) = idx(i+1). = index in cur_probes
            cur_pclass.boundaries.end_idx = cur_pclass.idx[which(cur_pclass.difference != 
                1)]  #
            cur_pclass.boundaries.start_idx = c(head(cur_pclass.idx, 
                1), cur_pclass.boundaries.start_idx)  #first region starts at the first idx in cur_probes with this class
            cur_pclass.boundaries.end_idx = c(cur_pclass.boundaries.end_idx, 
                tail(cur_pclass.idx, 1))  #last regeion ends at the last idx in cur_probes with this class
            
            cur_pclass.boundaries = data.frame(cur_pclass.boundaries.start_idx, 
                cur_pclass.boundaries.end_idx, cur_pclass.boundaries.end_idx - 
                  cur_pclass.boundaries.start_idx + 1)  #create data frame for easy access and add column with number of probes per region
            colnames(cur_pclass.boundaries) = c("start.idx", "end.idx", 
                "n_probes")
            
            # just ignore regions with too little probes
            n1 = which(cur_pclass.boundaries$n_probes < n_min)
            if (length(n1) > 0) {
                cur_pclass.boundaries = cur_pclass.boundaries[-n1, 
                  ]
            }
            
            if (length(cur_pclass.boundaries$n_probes) > 0) {
                # check with these regions if probes are very (>d) far appart; if
                # so, split and re-assess probe density
                for (cur_boundary in 1:dim(cur_pclass.boundaries)[1]) {
                  cur_boundary.idx = as.matrix(cur_pclass.boundaries$start.idx[cur_boundary]:cur_pclass.boundaries$end.idx[cur_boundary])  #take adjecent probes from current boundary (contains index in cur_probes)
                  cur_boundary.distance = as.matrix(diff(cur_probes[cur_boundary.idx, 
                    ]$position))  #calculate difference in positions between adjecent probes within current boundary
                  
                  # if distance between probes >= d --> split region
                  if (length(which(cur_boundary.distance >= d) > 
                    0)) {
                    cur_boundary.boundaries.start_idx = cur_boundary.idx[which(cur_boundary.distance >= 
                      d) + 1]
                    cur_boundary.boundaries.end_idx = cur_boundary.idx[which(cur_boundary.distance >= 
                      d)]
                    cur_boundary.boundaries.distance.start_idx = c(head(cur_boundary.idx, 
                      1), cur_boundary.boundaries.start_idx)
                    cur_boundary.boundaries.distance.end_idx = c(cur_boundary.boundaries.end_idx, 
                      tail(cur_boundary.idx, 1))
                  } else {
                    # OK, so no probes that are too far away from eachother --> copy
                    # data from cur_pclass.boundaries (=regions based on only
                    # annotation)
                    cur_boundary.boundaries.start_idx = cur_pclass.boundaries$start.idx[cur_boundary]
                    cur_boundary.boundaries.end_idx = cur_pclass.boundaries$end.idx[cur_boundary]
                    cur_boundary.boundaries.distance.start_idx = c(head(cur_boundary.idx, 
                      1))
                    cur_boundary.boundaries.distance.end_idx = c(tail(cur_boundary.idx, 
                      1))
                  }
                  
                  # just ignore regions with too little probes
                  np = cur_boundary.boundaries.distance.end_idx - 
                    cur_boundary.boundaries.distance.start_idx + 
                    1
                  n1 = which(np < n_min)
                  if (length(n1) > 0) {
                    # circumvent R behavior: if n1=empty--> all rows are discarded
                    # when using [-n1] which is exactly the oposite of what we want
                    cur_boundary.boundaries.distance.start_idx = cur_boundary.boundaries.distance.start_idx[-n1]
                    cur_boundary.boundaries.distance.end_idx = cur_boundary.boundaries.distance.end_idx[-n1]
                  }
                  
                  if (length(cur_boundary.boundaries.distance.start_idx) > 
                    0) {
                    lb = 1 + regionID
                    ub = length(cur_boundary.boundaries.distance.end_idx) + 
                      regionID
                    cur_boundary.boundaries = data.frame(cur_probes[cur_boundary.boundaries.distance.start_idx, 
                      ]$chr, cur_probes[cur_boundary.boundaries.distance.start_idx, 
                      ]$position, cur_probes[cur_boundary.boundaries.distance.end_idx, 
                      ]$position, cur_probes[cur_boundary.boundaries.distance.end_idx, 
                      ]$position - cur_probes[cur_boundary.boundaries.distance.start_idx, 
                      ]$position + 1, cur_boundary.boundaries.distance.start_idx, 
                      cur_boundary.boundaries.distance.end_idx, cur_boundary.boundaries.distance.end_idx - 
                        cur_boundary.boundaries.distance.start_idx + 
                        1, lb:ub, rep(cur_pclass, length(lb:ub)), 
                      rep(u_pclass[cur_pclass], length(lb:ub)))
                    
                    colnames(cur_boundary.boundaries) = c("chr", 
                      "start_bp", "end_bp", "length_bp", "start.idx", 
                      "end.idx", "n_probes", "regionID", "classID", 
                      "class")
                    
                    regionID = dim(cur_boundary.boundaries)[1] + 
                      regionID
                    regions = rbind(regions, cur_boundary.boundaries)  #merge
                    for (i in 1:dim(cur_boundary.boundaries)[1]) {
                      regions2probes[cur_probes$rowID[cur_boundary.boundaries$start.idx[i]:cur_boundary.boundaries$end.idx[i]], 
                        cur_boundary.boundaries$classID[i]] = cur_boundary.boundaries$regionID[i]
                    }
                  }
                }
            }
        }
    }
    colnames(regions2probes) = u_pclass
    
    # find exactly overlapping regions from different classes and
    # merge.
    regions_unique = unique(regions[, c(1, 2, 3, 4, 7)])
    k = dim(regions_unique)[1]
    
    regions_unique$regionID = ""
    regions_unique$regionIDall = ""
    regions_unique$ClassAll = ""
    
    for (i in 1:k) {
        cur_IDs = which(regions[, 1] == regions_unique[i, 1] & regions[, 
            2] == regions_unique[i, 2] & regions[, 3] == regions_unique[i, 
            3] & regions[, 4] == regions_unique[i, 4] & regions[, 
            7] == regions_unique[i, 5])
        cur_rows = regions[cur_IDs, ]
        regions_unique$regionID[i] = cur_rows$regionID[1]
        regions_unique$regionIDall[i] = paste(cur_rows$regionID, 
            collapse = ";")
        regions_unique$ClassAll[i] = paste(cur_rows$class, collapse = ";")
    }
    
    rownames(regions2probes) = rownames(m.v)
    rownames(probes) = rownames(m.v)
    return(list(boundaries = regions_unique, perprobe = regions2probes, 
        valid.probes = probes, valid.m = m.v, valid.beta = beta.v, 
        gs = gs))
}

testregion <- function(x, probes, m, b, n, dMth, do.format = FALSE)
{
    probe_rows = which(probes == x, arr.ind = TRUE)  #rows of probes in current regions
    probe_rows = probe_rows[, 1]
    out = calc_stats(probe_rows = probe_rows, probes, m, b, n, dMth, do.format = do.format)
    out
}

calc_stats <- function(probe_rows, probes, m, b, n, dMth, do.format = FALSE)
{
    tb = b[probe_rows, ]  #m and beta values in current region
    tm = m[probe_rows, ]
    np = dim(tb)[1]  #number of probes in current region
    
    i = combn(n, 2)
    mdb = apply(as.data.frame(tb[, i[1, ]] - tb[, i[2, ]]), 2, median)  #create 2 matrices with columns representing all possible pairwise comparisons and calculate median differences
    mdm = apply(as.data.frame(tm[, i[1, ]] - tm[, i[2, ]]), 2, median)
    max.mdm = max(abs(mdm))
    
    p = NA
    pairwise = rep(NA, dim(i)[2])
    
    if (max.mdm > dMth) {
        # if max median dm is large enough, test difference formally
        if (n == 2) {
            p = wilcox.test(tm[, 1], tm[, 2], paired = FALSE)$p.value
        } else if (n > 2) {
            g = rep(1:n, np)
            g = as.data.frame(g[order(g)])
            colnames(g) = "sample"
            g$m = c(apply(t(tm), 1, rbind))
            p = kruskal.test(m ~ sample, data = g)$p.value
            if (p <= 0.05) 
                {
                  for (j in 1:dim(i)[2]) {
                    pairwise[j] = wilcox.test(tm[, i[1, j]], tm[, 
                      i[2, j]], paired = FALSE)$p.value
                  }
                }  #If KW is sign (uncorrected), perform pairwise testing.
        }
    }
    out = c(apply(tb, 2, median), apply(tm, 2, median), mdb, mdm, 
        pairwise, max.mdm, p)
    
    if (do.format == TRUE) {
        out = t(out)
        cn1 = paste("beta.median", colnames(b), sep = ".")
        cn2 = paste("m.median", colnames(m), sep = ".")
        cn3 = paste("median.delta.beta", colnames(m)[i[1, ]], "minus", 
            colnames(m)[i[2, ]], sep = ".")
        cn4 = paste("median.delta.m", colnames(m)[i[1, ]], "minus", 
            colnames(m)[i[2, ]], sep = ".")
        cn5 = paste("pairwise.p", colnames(m)[i[1, ]], "vs", colnames(m)[i[2, 
            ]], sep = ".")
        colnames(out) = c(cn1, cn2, cn3, cn4, cn5, "max.abs.median.delta.m", 
            "p.value")
        out = as.data.frame(out)
        out$n.probes = length(probe_rows)
        out = t(out)
    }
    out
}

mod <- function(x, m)
{
    t1 <- floor(x/m)
    return(x - t1 * m)
}

matchsymbol <- function(l, str, sep = ";")
{
    # find complete matches of str in l, accepting that each row in l
    # has multiple entries separated by sep.
    l = as.character(l)
    match.exact = which(l == str)  #returns index in l
    pm = regexpr(str, l)  # returns > 1 if match
    match.partial = which(pm > 0)  #find probes which have at least the current class in their 
    match.partial.only = setdiff(match.partial, match.exact)  #exclude exact matches
    match.partial.only.exact = match.partial.only  #preprocess for loop
    
    if (!(length(match.partial.only) == 0)) {
        for (j in 1:length(match.partial.only)) {
            split = data.frame(strsplit(l[match.partial.only[j]], 
                sep))  #split partial match at sep and interpret the resulting strings separately
            if (length(which(split == str)) == 0) {
                # if none of the splitted strings matches str exactly, remove the
                # index from l from match.partial.only.exact
                match.partial.only.exact = match.partial.only.exact[-which(match.partial.only.exact == 
                  match.partial.only[j])]
            }
        }
    }
    match = c(match.exact, match.partial.only.exact) #only return rows of which at least 1 substring exactly matched str.
}

tune_parameters <- function(parameters, classes_gene, classes_island, 
    targetID, chr, position, m.v, beta.v, recode = 1, sep = ";", 
    gs, do.parallel = 0)
{
    message(paste("min_distance=", parameters[1], ", min_n=", parameters[2], sep = ""))
    # merge gene/island classes
    classes <- merge_classes(classes_gene, classes_island, recode = recode, sep = ";")
    
    message("Calculating the number of regions and associated probes for the requested set of parameters")
    if (do.parallel == FALSE) {
        results = apply(parameters, 1, function(parameters) tune_parameters_calc(parameters, 
            targetID, classes, chr, position, m.v, beta.v, recode = 1, 
            sep = ";", gs))
    }
    
    if (do.parallel == -1) {
        do.parallel = detectCores()
    }
    if (do.parallel > 0) {
        cl <- makeCluster(getOption("cl.cores", do.parallel))
        clusterExport(cl, c("tune_parameters_calc", "regionfinder", 
            "targetID", "classes", "chr", "position", "m.v", "beta.v", 
            "recode", "sep", "gs"), envir = environment())
        results = parApply(cl, parameters, 1, function(parameters) tune_parameters_calc(parameters, 
            targetID, classes, chr, position, m.v, beta.v, recode, 
            sep, gs))
        stopCluster(cl)
    }
    
    results = t(results)
    
    colnames(results) = c("min_distance", "min_n", "n.regions", "n.valid.probes", 
        "n.probes.included")
    return(results)
}

tune_parameters_calc <- function(parameters, targetID, classes, chr, 
    position, m.v, beta.v, recode = 1, sep = ";", gs)
{
    # find potential regions of interest
    regions <- regionfinder(targetID = targetID, chr = chr, position = position, 
        classes$pclass_recoded, classes$no.pclass, classes$u_pclass, 
        d = parameters[1], m.v = m.v, beta.v = beta.v, n_min = parameters[2], 
        debug = FALSE, gs = gs)
    # summarize the results
    results = matrix()
    results[1] = parameters[1]
    colnames(results) = "d"
    results[2] = parameters[2]
    results[3] = dim(regions$boundaries)[1]
    results[4] = dim(regions$valid.probes)[1]
    tmp = regions$perprobe
    j = which(is.na(regions$perprobe))
    tmp[j] = 0
    results[5] = length(which(rowSums(tmp) > 0))
    return(results)
}

export_data <- function(tested, regions, th = 0.05, annotate.relevant = FALSE, 
    annotate.significant = TRUE, FigsNotRelevant = FALSE, min_n = 4, 
    min_dM = 1.4, min_distance = 200, margin = 10000, clr = NA, method = "fdr", 
    experiment.name, debug = FALSE)
{
    dir.create(file.path(getwd(), experiment.name))
    
    dir.create(file.path(paste(getwd(), experiment.name, sep = "/"), "figures"))
    old.wd = getwd()
    
    setwd(file.path(getwd(), experiment.name))
    path = "figures"
    message("Preparing export")
    # add columns to tested with links to the regions in UCSC and
    # ENSEMBL
    linkEnsembl_pre = "<a href=\"http://www.ensembl.org/Homo_sapiens/Location/View?r="
    linkEnsembl_post = "\" target=\"_blank\">ENSEMBL</a>"
    tested$LinkEnsembl = paste(linkEnsembl_pre, tested$chr, ":", tested$start_bp, "-", tested$end_bp, linkEnsembl_post, sep = "")
    
    linkUCSC_pre = "<a href=\"http://genome-euro.ucsc.edu/cgi-bin/hgTracks?position=chr"
    linkUCSC_post = "&hgsid=192199020&knownGene=pack&hgFind.matches=uc004dqr.3,\" target=\"_blank\">UCSC</a>"
    tested$LinkUCSC = paste(linkUCSC_pre, tested$chr, ":", tested$start_bp, "-", tested$end_bp, linkUCSC_post, sep = "")
    
    # export figures / pdfs and add annotation to tested
    tested$GeneSymbols_exact = "Not queried"
    tested$GeneSymbols_margin = "Not queried"
    tested$Figure = NA
    tested$Statistics = NA
    
    # regions not of interest (not significant or not relevant
    # (median dM < threshold))
    if (FigsNotRelevant == TRUE) {
        message("Generating images for not relevant regions")
        r = tested[which(is.na(tested$p.value.adjusted)), "regionID"]
        for (i in 1:length(r)) {
            annot = plot_annotate_region(tested, regions, margin = margin, 
                regionID = r[i], clr = clr, annotate = FALSE, scores = FALSE, 
                path = path)
            cr = which(tested$regionID == r[i])
            tested$Figure[cr] = paste("<a href=\"./", path, "/", 
                tested$regionID[cr], ".pdf\" target=\"_blank\">PDF</a>", 
                sep = "")
            tested$Statistics[cr] = paste("<a href=\"./", path, "/", 
                tested$regionID[cr], ".tsv\" target=\"_blank\">STATS</a>", 
                sep = "")
            write.table(file = paste("./", path, "/", tested$regionID[cr], 
                ".tsv", sep = ""), t(tested[cr, ]), sep = "\t", col.names = FALSE, 
                row.names = TRUE)
        }
    }
    
    # relevant, but not significant regions
    r = tested[which(!is.na(tested$p.value.adjusted) & tested$p.value.adjusted > th), "regionID"]
    if (length(r) > 0) {
        message("Generating images for relevant regions (if annotation is requested, this can take quite long).")
        for (i in 1:length(r)) {
            annot = plot_annotate_region(tested, regions, margin = margin, 
                regionID = r[i], clr = clr, annotate = annotate.relevant, 
                scores = FALSE, path = path)
            cr = which(tested$regionID == r[i])
            tested$Figure[cr] = paste("<a href=\"./", path, "/", 
                tested$regionID[cr], ".pdf\" target=\"_blank\">PDF</a>", sep = "")
            tested$Statistics[cr] = paste("<a href=\"./", path, "/", 
                tested$regionID[cr], ".tsv\" target=\"_blank\">STATS</a>", sep = "")
            if (annotate.relevant == TRUE) {
                tested$GeneSymbols_exact[cr] = annot$symbols_exact
                tested$GeneSymbols_margin[cr] = annot$symbols_margin
            }
            write.table(file = paste("./", path, "/", tested$regionID[cr], 
                ".tsv", sep = ""), t(tested[cr, ]), sep = "\t", col.names = FALSE, 
                row.names = TRUE)
        }
    }
    
    # regions interest (significant & relevant (median dM >=
    # threshold))
    r = tested[which(tested$p.value.adjusted <= th), "regionID"]
    if (length(r) > 0) {
        message("Generating images for significant regions (if annotation is requested, this can take quite long).")
        for (i in 1:length(r)) {
            annot = plot_annotate_region(tested, regions, margin = margin, 
                regionID = r[i], clr = clr, annotate = annotate.significant, 
                scores = FALSE, path = path)
            cr = which(tested$regionID == r[i])
            if (annotate.significant == TRUE) {
                tested$GeneSymbols_exact[cr] = annot$symbols_exact
                tested$GeneSymbols_margin[cr] = annot$symbols_margin
            }
            tested$Figure[cr] = paste("<a href=\"./", path, "/", 
                tested$regionID[cr], ".pdf\" target=\"_blank\">PDF</a>", sep = "")
            tested$Statistics[cr] = paste("<a href=\"./", path, "/", 
                tested$regionID[cr], ".tsv\" target=\"_blank\">STATS</a>", sep = "")
            write.table(file = paste("./", path, "/", tested$regionID[cr], 
                ".tsv", sep = ""), t(tested[cr, ]), sep = "\t", col.names = FALSE, row.names = TRUE)
        }
    }
    
    message("Exporting tables (HTML & CSV)")
    # for easy sorting and selection...
    tested$links = paste(tested$Figure, tested$Statistics, tested$LinkEnsembl, tested$LinkUCSC, sep = "</br>")
    tested$links = gsub("NA</br>", "", tested$links)  #remove empty entries
    
    # merge some columns to limit the number of columns in the output
    tested$gene.symbols = paste(tested$GeneSymbols_exact, " (margin: ", tested$GeneSymbols_margin, ")", sep = "")
    tested$gene.symbols[which(tested$gene.symbols == "Not queried (margin: Not queried)")] = "Not queried"
    
    # format data frame for output to html & tsv
    i = which(regexpr("beta.median", colnames(tested)) > 0)  #columns with median beta values
    tested_4html = tested[, c("chr", "start_bp", "end_bp", "links", 
        "length_bp", "n_probes", "regionID", "ClassAll", colnames(tested)[i], 
        "gene.symbols", "max.abs.median.delta", "p.value", "p.value.adjusted")]
    colnames(tested_4html)[1:8] = c("Chr", "Start", "End", "Links", "Length", "n", "ID", "Class")
    colnames(tested_4html) = gsub("beta.median.", "", colnames(tested_4html))
    colnames(tested_4html)[(dim(tested_4html)[2] - 3):dim(tested_4html)[2]] = c("Gene.Symbol", "dM", "p", "p.adj")
    tested_4html = tested_4html[order(tested_4html$p.adj, -tested_4html$dM), ]
    options(scipen = 4)
    options(digits = 2)
    
    leg = paste("Chromosomal positions are indicated in bp. N indicates the number of probes in a region. ID indicates the region ID. Values in the columns per sample indicate median beta values. dM indicates the largest median difference (absolute) between any of the sample pairs. p indicates uncorrected p-value from Mann-Whitney U test (n=2) or Kruskal Wallis test (n>2). p.adj denotes the multiple testing corrected p-value (method:", 
        method, "). Gene symbols of overlapping transcripts are listed for the exact region and within a margin of ", 
        margin, " bp of the region.", sep = "")
    caption = paste("DMRforPairs output generated on ", date(), ". Identified regions were set to contain at least ", 
        min_n, " probes with a maximum distance of ", min_distance, 
        " bp between individual probes (n=", length(tested[, 1]), 
        "). Regions in which median methylation levels (M-values) between the samples differed at least |", 
        min_dM, "| (n=", length(which(tested$p.value.adjusted <= 
            1)), ") (=relevant) were tested for statistical significance (significant: p<", 
        th, "; multiple testing adjusted, n=", length(which(tested$p.value.adjusted <= 
            th)), "). The folowing samples were studied: ", paste(colnames(regions$valid.m), 
            collapse = ","), ". ", leg, sep = "")
    
    # export all regions (no thumbnails)
    HTML(tested_4html, file = "all.html", row.names = FALSE, align = "left", 
        Border = NULL, innerBorder = 1, append = FALSE, caption = caption, 
        captionalign = "top", digits = 3, big.mark = "", big.interval = 3, 
        decimal.mark = ".")
    write.table(tested_4html, "all.tsv", sep = "\t", row.names = FALSE)
    
    # export relevant regions (dM sufficiently large)
    tested_4html_selected = tested_4html[which(tested_4html$p.adj <= 1), ]
    if (dim(tested_4html_selected)[1] > 0) {
        tested_4html_selected = tested_4html_selected[order(-tested_4html_selected$dM), ]
        tested_4html_selected = cbind("", tested_4html_selected)
        colnames(tested_4html_selected)[1] = "Thumbnail"  #add thumbnail column
        tested_4html_selected$Thumbnail = paste("<img src=\"./", path, "/", tested_4html_selected$ID, ".png\" height=\"63\" width=\"125\">", sep = "")
        HTML(tested_4html_selected, file = "relevant.html", row.names = FALSE, 
            align = "left", Border = NULL, innerBorder = 1, append = FALSE, 
            caption = caption, captionalign = "top", digits = 3, 
            big.mark = "", big.interval = 3, decimal.mark = ".")
        write.table(tested_4html_selected, "relevant.tsv", sep = "\t", row.names = FALSE)
    }
    # export significant DMRs
    tested_4html_selected = tested_4html[which(tested_4html$p.adj <= th), ]
    if (dim(tested_4html_selected)[1] > 0) {
        tested_4html_selected = tested_4html_selected[order(-tested_4html_selected$dM), ]
        tested_4html_selected = cbind("", tested_4html_selected)
        colnames(tested_4html_selected)[1] = "Thumbnail"  #add thumbnail column
        tested_4html_selected$Thumbnail = paste("<img src=\"./", path, "/", tested_4html_selected$ID, ".png\" height=\"63\" width=\"125\">", sep = "")
        HTML(tested_4html_selected, file = "significant.html", row.names = FALSE, 
            align = "left", Border = NULL, innerBorder = 1, append = FALSE, 
            caption = caption, captionalign = "top", digits = 3, 
            big.mark = "", big.interval = 3, decimal.mark = ".")
        write.table(tested_4html_selected, "significant.tsv", sep = "\t", row.names = FALSE)
    }
    
    setwd(old.wd)
    return(tested)
}

plot_annotate_region <- function(tested, regions, margin = 10000, 
    regionID, clr = NA, annotate = TRUE, scores = TRUE, path)
{
    # wrapper for plot_annotate_probes to plot / find annotation for
    # one region
    path = paste("./", path, "/", sep = "")
    cur_row = which(tested$regionID == regionID)
    title_x = paste("RegionID: ", regionID, ", chr", tested$chr[cur_row], 
        ":", tested$start_bp[cur_row], "-", tested$end_bp[cur_row], sep = "")
    probe_rows <- which(regions$perprobe == regionID, arr.ind = TRUE)
    probe_rows <- probe_rows[, 1]
    annot <- plot_annotate_probes(regions = regions, title_x = title_x, 
        probe_rows = probe_rows, margin = margin, ID = regionID, 
        clr = clr, annotate = annotate, scores = scores, path = path)
    return(annot)
}

plot_annotate_gene <- function(gs, regions, margin = 10000, ID, clr = NA, 
    annotate = TRUE, path)
{
    # wrapper for plot_annotate_probes to plot / find annotation for
    # one gene (gene symbol)
    path = paste("./", path, "/", sep = "")
    probe_rows = matchsymbol(regions$gs, gs, sep = ";")
    
    if (length(probe_rows) == 0) {
        message(paste("Gene Symbol ", gs, " not found", sep = ""))
        annot = matrix()
    } else {
        chr = regions$valid.probes$chr[probe_rows[1]]
        st = min(regions$valid.probes$position[probe_rows])
        ed = max(regions$valid.probes$position[probe_rows])
        title_x = paste(ID, ", chr", chr, ":", st, "-", ed, sep = "")
        
        annot = plot_annotate_probes(regions = regions, title_x = title_x, 
            probe_rows = probe_rows, margin = margin, ID = ID, clr = clr, 
            annotate = annotate, scores = TRUE, path = path)
    }
    return(annot)
}

plot_annotate_custom_region <- function(chr, st, ed, regions, margin = 10000, 
    ID = "CustomRegion", clr = NA, annotate = TRUE, path)
{
    path = paste("./", path, "/", sep = "")
    # wrapper for plot_annotate_probes to plot / find annotation for
    # a custom genomic region
    probe_rows = which(regions$valid.probes$chr == chr & regions$valid.probes$position <= 
        ed & regions$valid.probes$position >= st)
    title_x = paste(ID, ", chr", chr, ":", st, "-", ed, sep = "")
    if (length(probe_rows) == 0) {
        message(paste("No valid probes present in the specified region", 
            sep = ""))
        annot = matrix()
    } else {
        annot = plot_annotate_probes(regions = regions, title_x = title_x, 
            probe_rows = probe_rows, margin = margin, ID = ID, clr = clr, 
            annotate = annotate, scores = TRUE, path = path)
    }
    return(annot)
}

plot_annotate_probes <- function(regions, title_x, probe_rows, margin = 10000, 
    ID = NA, clr = NA, annotate = TRUE, scores = NA, path) # plot_annotate_probes
{
    # Little quirk of Gviz that needs to be compensated for to match
    # the colors in all graphs...
    regions$valid.m = regions$valid.m[, order(colnames(regions$valid.m))]
    regions$valid.beta = regions$valid.beta[, order(colnames(regions$valid.beta))]
    
    # plot methylation values and genomic annotation info and output
    # gene symbols + gene region within and near probes (region) of
    # interest
    ns = dim(regions$valid.beta)[2]
    chro = regions$valid.probes$chr[probe_rows[1]]
    st = min(regions$valid.probes$position[probe_rows])
    ed = max(regions$valid.probes$position[probe_rows])
    
    # If custom region is requested, deliver statistics as well.
    if (scores == TRUE) {
        probes <- regions$perprobe
        m <- regions$valid.m
        b <- regions$valid.beta
        n <- ns
        dMth <- 0  #this is not happening when regions are selected by DMRforPairs, but in cases where custom regions are requested by the user: always return test values then.
        scores = calc_stats(probe_rows, probes, m, b, n, dMth, do.format = TRUE)
    }
    
    if (!length(clr) == ns) {
            clr = rainbow(ns)
        }  #Did the user specify colors? If not, pick some from the rainbow pallet.
    
    # mini fig for use in html overview
    pos = regions$valid.probes$position[probe_rows]
    png(paste(path, ID, ".png", sep = ""), width = 500, height = 250)
    # plot beta values
    par(mai = c(0, 0, 0, 0))
    plot(sort(rep(pos, ns)), t(regions$valid.beta[probe_rows[order(pos)], 
        ]), col = clr, lty = 2.5, cex = 2.5, pch = c(19), ylim = c(-0.1, 
        1.1), xaxt = "n", yaxt = "n", ann = FALSE, frame = FALSE)
    lines(extendrange(regions$valid.probes$position[probe_rows]), c(0, 0), lty = "dashed")
    lines(extendrange(regions$valid.probes$position[probe_rows]), c(0.25, 0.25), lty = "dashed")
    lines(extendrange(regions$valid.probes$position[probe_rows]), c(0.5, 0.5), lty = "dashed")
    lines(extendrange(regions$valid.probes$position[probe_rows]), c(0.75, 0.75), lty = "dashed")
    lines(extendrange(regions$valid.probes$position[probe_rows]), c(1, 1), lty = "dashed")
    dev.off()
    
    # pdf with detailed plots and annotation
    pdf(paste(path, ID, ".pdf", sep = ""), width = 10, height = 10)
    
    par(mfrow = c(2, 1))
    par(xpd = TRUE)
    
    # plot M-values
    up = max(regions$valid.m[probe_rows, ]) + 0.5
    down = min(regions$valid.m[probe_rows, ]) - 0.5
    
    i = combn(ns, 2)
    
    pos = regions$valid.probes$position[probe_rows]
    
    for (j in 0:dim(i)[2]) {
        if (j == 0) {
            # overall plot
            m = regions$valid.m[probe_rows[order(pos)], ]
            b = regions$valid.beta[probe_rows[order(pos)], ]
            cur_clr = clr
            s.names = colnames(regions$valid.beta)
            k = ns
            go = TRUE
        } else if (ns > 2) {
            # pairwise plots if > 2 samples
            m = regions$valid.m[probe_rows, c(i[, j])]
            m = m[order(pos), ]
            b = regions$valid.beta[probe_rows, c(i[, j])]
            b = b[order(pos), ]
            cur_clr = clr[c(i[, j])]
            s.names = colnames(regions$valid.beta)[c(i[, j])]
            k = 2
            go = TRUE
        } else if (ns == 2) 
            {
                go = FALSE
            }  #prevent pairwise plots if n=2
        
        if (go == TRUE) {
            plot(sort(rep(pos, k)), t(m), col = cur_clr, lwd = 4, 
                cex = 4, pch = c(21), main = paste(title_x, "-M_values", 
                  sep = ""), xlab = "position (bp)", ylab = "M", 
                ylim = c(down, up), xlim = c(st - 0.17 * (ed - st), ed))
            legend("topleft", legend = s.names, col = cur_clr, pch = rep(21, ns), cex = 1)  #,inset=c(0,-0.2)
            
            # plot beta values
            plot(sort(rep(pos, k)), t(b), col = cur_clr, lwd = 4, 
                cex = 4, pch = c(21), main = paste(title_x, "-Beta_values", 
                  sep = ""), xlab = "position (bp)", ylab = "Beta", 
                ylim = c(0, 1), xlim = c(st - 0.17 * (ed - st), ed), 
                axes = FALSE)
            axis(2, c(0, 0.25, 0.5, 0.75, 1), labels = TRUE)
            axis(1, labels = TRUE)
            axis(4, seq(0, 1, 0.05), labels = FALSE)
            box(which = "plot", lty = "solid")
            legend("topleft", legend = s.names, col = cur_clr, pch = rep(21, ns), cex = 1)  #,inset=c(0,-0.2)
            
            # points(regions$valid.probes$position[probe_rows],regions$valid.beta[probe_rows,1],col=clr[2],lwd=6,cex=4,pch
            # = c(21))
            lines(extendrange(regions$valid.probes$position[probe_rows]), c(0, 0), lty = "dashed")
            lines(extendrange(regions$valid.probes$position[probe_rows]), c(0.25, 0.25), lty = "dashed")
            lines(extendrange(regions$valid.probes$position[probe_rows]), c(0.5, 0.5), lty = "dashed")
            lines(extendrange(regions$valid.probes$position[probe_rows]), c(0.75, 0.75), lty = "dashed")
            lines(extendrange(regions$valid.probes$position[probe_rows]), c(1, 1), lty = "dashed")
        }
    }
    
    if (annotate == TRUE) {
        options(ucscChromosomeNames = FALSE)
        # annotation track (exact region)
        
        bmTrack <- BiomartGeneRegionTrack(start = st, end = ed, chromosome = chro, 
            genome = "hg19", showId = TRUE, background.title = "white", 
            name = "", fontsize = 20, col = "black")
        
        symbols_exact = paste(unique(symbol(bmTrack)), collapse = "; ")
        rm(bmTrack)
        
        # annotation track (region+margin)
        prm = which(regions$valid.probes$chr == chro & regions$valid.probes$position <= 
            ed + margin & regions$valid.probes$position >= st - margin)  #find all probes in region +-margin
        
        
        bmTrack <- BiomartGeneRegionTrack(start = st - margin, end = ed + 
            margin, chromosome = chro, genome = "hg19", showId = TRUE, 
            background.title = "white", name = "", fontsize = 20, 
            col = "black", collapseTranscripts = FALSE)  #,col.line = NULL, stackHeight = 0.3,fontsize=11, #collapseTranscripts=TRUE    
        symbols_margin = paste(unique(symbol(bmTrack)), collapse = "; ")
        
        # Convert al tracks to 'chrX' like chr identifiers.
        # BiomartGeneRegionTrack based on e! which is using 1, 2, 3, ...
        # X IdeogramTrack based on UCSC which is using chr1, chr2, chr3,
        # ... chrX
        seqlevels(ranges(bmTrack)) <- sprintf("chr%s", seqlevels(ranges(bmTrack)))
        
        chromosome(bmTrack) = sprintf("chr%s", chro)
        
        # ideogram & genome axis
        itrack <- IdeogramTrack(genome = "hg19", chromosome = sprintf("chr%s", 
            chro), fontsize = 20)
        gtrack <- GenomeAxisTrack(fontsize = 20)
        
        # AnnotationTrack to indicate ROI
        st <- c(st)
        ed <- c(ed)
        strand <- c("*")
        gr <- c("ROI")
        
        annTrack <- AnnotationTrack(start = st, end = ed, strand = strand, 
            chromosome = sprintf("chr%s", chro), genome = "hg19", 
            feature = "ROI", group = gr, id = paste(ID), name = "generic annotation", 
            stacking = "squish", background.title = "white", col = "black", 
            fill = "black", showFeatureId = FALSE, fontcolor = "black")
        
        for (j in 0:dim(i)[2]) {
            if (j == 0) {
                # overall plot
                genes <- GRanges(seqnames = Rle(c(paste("chr", regions$valid.probes$chr[prm], 
                  sep = ""))), ranges = IRanges(start = regions$valid.probes$position[prm], 
                  end = regions$valid.probes$position[prm], names = regions$valid.probes$TagetID), 
                  strand = Rle(rep("*", length(prm))), regions$valid.beta[prm, 
                    ])
                cur_clr = clr
                s.names = colnames(regions$valid.beta)
                k = ns
                go = TRUE
            } else if (ns > 2) {
                # pairwise plots if > 2 samples
                genes <- GRanges(seqnames = Rle(c(paste("chr", regions$valid.probes$chr[prm], 
                  sep = ""))), ranges = IRanges(start = regions$valid.probes$position[prm], 
                  end = regions$valid.probes$position[prm], names = regions$valid.probes$TagetID), 
                  strand = Rle(rep("*", length(prm))), regions$valid.beta[prm, 
                    c(i[, j])])
                cur_clr = clr[c(i[, j])]
                s.names = colnames(regions$valid.beta)[c(i[, j])]
                k = 2
                go = TRUE
            } else if (ns == 2) 
                {
                  go = FALSE
                }  #prevent pairwise plots if n=2
            
            if (go == TRUE) {
                dTrack <- DataTrack(genes, name = "\U03B2", groups = s.names, 
                  legend = TRUE, , col = cur_clr, cex = 2, pch = c(19), 
                  fill = "transparent", fontsize = 16, ylim = c(0, 
                    1.001), fontsize.legend = 20, fontsize.title = 16, 
                  col.axis = "black", col.title = "black", background.title = "white", 
                  col.grid = "black", col.name = "black", type = c("p", 
                    "g"), col.grid = "grey75", v = 0, lwd = 4)
                
                plotTracks(list(itrack, dTrack, annTrack, gtrack, 
                  bmTrack), chromosome = sprintf("chr%s", chro), 
                  from = st - margin, to = ed + margin, extend.left = 0)
            }
        }
        
        # clean up
        rm(annTrack)
        rm(bmTrack)
        rm(itrack)
        rm(gtrack)
    } else {
        symbols_exact = ""
        symbols_margin = ""
    }
    dev.off()
    return(list(symbols_exact = symbols_exact, symbols_margin = symbols_margin, 
        scores = scores))
}
