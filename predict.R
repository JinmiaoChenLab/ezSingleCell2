args.input=commandArgs(T)

workdir = args.input[1]
refDataset = args.input[2]
library(pbmcapply)
#library(Rfast)
library(stringr)
library(RJSONIO)
#library(Seurat)


predict.ct = function(input.data, ref.data, ncore = 1){

  genes = intersect(rownames(ref.data), rownames(input.data))
  #print(genes)
  input.data = input.data[genes,]
  #print(input.data)
  ref.data = ref.data[genes,]
  #print(ref.data)

  ############## First round: use all gene ####################

  predicted.cell = pbmclapply(
    1:ncol(input.data),
    function(j) {
      predicted = as.numeric(
        apply(ref.data, 2, function(i) {
          cor(as.numeric(i), as.numeric(input.data[, j]), method = "spearman", use="complete.obs")
        })
      )
      return(predicted)
    }, mc.cores = ncore
  )

  #print (predicted.cell)

  predicted.cell = do.call(cbind, predicted.cell)
  rownames(predicted.cell) = colnames(ref.data)

  ct = apply(predicted.cell, 2, function(x) {
    return(as.numeric(which(x > 0.6 | rank(-x) <=5)))
  })

  #print (ct)
  predicted.cell = pbmclapply(
    1:ncol(input.data), function(i) {
      ref = ref.data[,ct[[i]]]
      print(ref)
      g = which(rank(-rowVars(as.matrix(ref))) <= 3000)
      print(g)
      ref = ref[g,]
      print(ref)
      input = input.data[g,i]
      predict = apply(ref, 2, function(i) {
        cor(as.numeric(i), input, method = "spearman", use="complete.obs")
      })
      return(c(names(sort(predict, decreasing = T)[1:3]), as.numeric(sort(predict, decreasing = T)[1:3])))
    }, mc.cores = ncore
  )

  #print (predicted.cell)


  predicted.cell = do.call(rbind, predicted.cell)

  cell.type.meta = fromJSON("https://www.immunesinglecell.org/api/cell_type/getReferenceCellType")
  cell.type.meta = do.call(rbind, cell.type.meta)
  cell.type.meta = data.frame(cell.type.meta)
  rownames(cell.type.meta) = paste0(cell.type.meta$dataset, ".", cell.type.meta$originalName)
  for (i in 1:10) {
    cell.type.meta[,i] = as.character(cell.type.meta[,i])
  }

  ontology = fromJSON("https://www.immunesinglecell.org/api/cell_type/getOntology")
  ontology = do.call(rbind, ontology)
  ontology = data.frame(ontology)
  ontology$idx = as.character(ontology$idx)
  ontology$parent = as.character(ontology$parent)
  ontology$children = as.character(ontology$children)

  res = data.frame(original_name = predicted.cell[,1], score = predicted.cell[,4], dataset = NA, standardized_name = NA,
                   publication = NA, tissue = NA, sample_status = NA, hierarchical = NA, cluster_number = 1:nrow(predicted.cell))

  for (i in 1:dim(res)[1]) {
    res[i, "dataset"] = cell.type.meta[res$original_name[i],"dataset"]
    res[i, "publication"] = cell.type.meta[res$original_name[i],"publication"]
    res[i, "tissue"] = cell.type.meta[res$original_name[i],"tissue"]
    res[i, "sample_status"] = cell.type.meta[res$original_name[i],"sampleStatus"]
    res[i, "standardized_name"] = cell.type.meta[res$original_name[i],"standardizedName"]

    ot = res[i, "standardized_name"]
    o = ontology$parent[which(ontology$children == ot)]
    while (o != "cell type") {
      ot = paste0(ot, ";", o)
      o = ontology$parent[which(ontology$children == o)]
    }

    res[i, "hierarchical"] = ot
    res[i, "original_name"] = cell.type.meta[res$original_name[i],"originalName"]
  }


  return(res)
}
