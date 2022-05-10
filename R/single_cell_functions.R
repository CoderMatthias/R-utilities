#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Custom Single Cell Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Load multiple 10x CellPlex output in Seurat object list
#'
#' This function takes information from a sample info dataframe containing the sample names, location of cellranger multi outputs, and metadata columns and loads each sample in as a Seurat object.
#' @param sample_info A sample information datafram containing; sample names, paths to cellranger multi outputs, and as many parameter columns as wanted to add to metadata
#' @param sample_name_column Name of column that contains the sample names
#' @param sample_location_column Name of column that contains the path to the cellranger multi outputs
#' @param metadata_columns Vector containing all the names of the columns that need to be added as meta data
#' @return Returns a list of Seurat objects to be used for downstream integration
#' @keywords Seurat
#' @export
#' @examples
#' plotPCAmulti()

SC_load_CellPlex <- function(sample_info, sample_name_column, sample_location_column, metadata_columns) {
  
  sample_list <- list()
  for (index in 1:nrow(sample_info)){  
    print(sample_info[[sample_name_column]][index])
    
    # Load 10x outputs
    sample.data <- Read10X(data.dir = as.character(sample_info[[sample_location_column]][index]))
    
    # check for duplicated barcodes(column names) within a dataset
    any(duplicated(x = colnames(x = sample.data))) %>% print()
    
    # Rename columns to avoid trouble when merging data sets (since identical columns between KO and WT)
    colnames(x = sample.data$`Gene Expression`) <- paste(as.character(sample_info[[sample_name_column]][index]), 
                                                         colnames(x = sample.data$`Gene Expression`), sep = '_')
    
    ### Initialize the Seurat object with the raw, non-normalized data
    # Keep all cells with at least 200 detected genes.
    sample <- CreateSeuratObject(counts = sample.data$`Gene Expression`, 
                                 min.cells = 0, 
                                 min.features = 200, 
                                 project = sample_info[[sample_name_column]][index])
    
    ### add treatment/group status to objects
    for(meta in metadata_columns) {
      sample[[meta]] <- sample_info[[meta]][index]
    }
    sample_list <- append(sample_list, sample)
  }
  return(sample_list)
}

