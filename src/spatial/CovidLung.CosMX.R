### the processed seurat object is too big to save it here, so please download it from our google drive through this link:
#  https://drive.google.com/file/d/10GO_GYqn2ZceVWcCNZS5WGxvuHIMWSxT/view?usp=drive_link

### setting and loading
library(Seurat)
library(future)
# plan("multisession", workers = 10)
library(ggplot2)
library(magrittr)
# install.packages("scCustomize")
library(scCustomize)
library(RColorBrewer)
library(writexl)
library(cli)
library(dplyr)
library(tidyr)
library(ggsignif)

### intermediate solution for FOV-optimized
if(T){
  #' @importFrom magrittr %>% %<>%
  NULL
  
  #' Intermediate solution to \code{subset()}:
  #' subset FOVs/centroids if selected cells are NOT found in each FOV
  #' NOTE: some code parts and args are taken from SeuratObject
  
  #' Function params/args:
  #' @param object An S4 object or A \code{FOV} object
  #' @param subset Logical expression indicating features/variables to keep
  #' @param cells A vector of cells to keep; if \code{NULL}, defaults to all cells
  #' @param idents A vector of identity classes to keep
  #' @param Update.slots If to update slots of an object
  #' @param Update.object If to update final object, default to TRUE.
  #' @param ... Arguments passed to \code{subset()} and other methods
  
  
  subset_opt <- function(
	object = NULL, 
	subset, 
	cells = NULL, 
	idents = NULL, 
	Update.slots = TRUE,
	Update.object = TRUE,
	...)
  {
	
	if (Update.slots) { 
	  message("Updating object slots..")
	  object %<>% UpdateSlots()
	}
	
	message("Cloing object..")
	obj_subset <- object
	
	# sanity check - use only cell ids (no indices)
	if (all(is.integer(cells))) { 
	  cells <- Cells(obj_subset)[cells]
	}
	
	if (!missing(subset) || !is.null(idents)) {
	  message("Extracting cells matched to `subset` and/or `idents`")
	}
	
	if (class(obj_subset) == "FOV") {
	  message("object class is `FOV` ")
	  cells <- Cells(obj_subset)
	} else if (!class(obj_subset) == "FOV" && !missing(subset)) {
	  subset <- enquo(arg = subset)
	  # cells to keep in the object
	  cells <-
		WhichCells(object = obj_subset, 
				   cells = cells,
				   idents = idents,
				   expression = subset,
				   return.null = TRUE, ...)
	} else if (!class(obj_subset) == "FOV" && !is.null(idents)) {
	  cells <-
		WhichCells(object = obj_subset, 
				   cells = cells,
				   idents = idents,
				   return.null = TRUE, ...)
	} else if (is.null(cells)) {
	  cells <- Cells(obj_subset)
	}
	
	# added support for object class `FOV`
	if (class(obj_subset) == "FOV") {
	  message("Matching cells for object class `FOV`..")
	  cells_check <- any(obj_subset %>% Cells %in% cells)
	} else { 
	  # check if cells are present in all FOV
	  message("Matching cells in FOVs..")
	  cells_check <-
		lapply(Images(obj_subset) %>% seq, 
			   function(i) { 
				 any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells) 
			   }) %>% unlist
	}
	
	if (all(cells_check)) { 
	  message("Cell subsets are found in all FOVs!", "\n",
			  "Subsetting object..")
	  obj_subset %<>% base::subset(cells = cells, idents = idents, ...)
	} else { 
	  # if cells are present only in one or several FOVs:
	  # subset FOVs
	  fovs <- 
		lapply(Images(obj_subset) %>% seq, function(i) {
		  if (any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)) {
			message("Cell subsets are found only in FOV: ", "\n", Images(obj_subset)[i])
			message("Subsetting Centroids..")
			base::subset(x = obj_subset[[Images(obj_subset)[i]]], cells = cells, idents = idents, ...)
		  }
		}) 
	  # replace subsetted FOVs, and remove FOVs with no matching cells
	  message("Removing FOVs where cells are NOT found: ", "\n", 
			  paste0(Images(object)[which(!cells_check == TRUE)], "\n"), "\n",
			  "Subsetting cells..")
	  for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }  
	  
	}
	
	# subset final object
	obj_subset %<>% base::subset(cells = cells, ...)
	
	if (Update.object && !class(obj_subset) == "FOV") { 
	  message("Updating object..")
	  obj_subset %<>% UpdateSeuratObject() }
	
	message("Object is ready!")
	return(obj_subset)
	
  }
}

### load raw data
obj <- readRDS("seuratObj.rds")
raw_data <- getwd()

### load markers
if(T){
  marker_dir <- file.path(raw_data,paste0("Marker.CosMX.xlsx"))
  sep <- ","
  library(readxl)
  ### merge the markers and plot violin and dot plot together
  library(readxl)
  all_markers <- read_excel(marker_dir,col_names = F)
  colnames(all_markers) <- c("cell_type","marker")
  markers.all <- c()
  for( i in 1:nrow(all_markers)){ ## i<-1
    cell_type <- all_markers[i,1]
    markers <- all_markers[i,2]
    markers <- unlist(strsplit(markers$marker, sep))
    markers <- gsub(" ","",markers)
    markers.all <- c(markers.all,markers) 
  }
  
}

col <- c("#66C5CC","#F89C74","#DCB0F2","#87C55F","#2F8AC4","#764E9F", "#DFFF00",
         "#0B775E","#9A8822","#7fcdbb","#dd1c77","#4D61C6","#0CD33F", "#696969" )

DimPlot(obj,reduction = "umap",label = T,label.size = 4,cols = col,label.color="black",pt.size = 1,alpha = 0.8,group.by = "cluV7")+
  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )+
  scale_color_manual( values = col,
                      labels = c("0:LC-1","1:Macro","2:Fibro","3:SMC","4:Endo","5:AT2","6:Mono","7:Plasma","8:RBC",
                                 "9:LC-2","10:AT1","11:Basal","12:Nutrophil","13:Unknown") )+
  theme( plot.title = element_blank())

DotPlot(obj, features = markers.all)

## dot plot for 1: macrophages and 6: monocytes using gene list of “pro-inflammatory associated genes” in excel
if(T){
  ## read pro-inflammatory associated genes
  marker_dir <- file.path(raw_data,"proinflammatory_related_genes.xlsx")
  library(readxl)
  all_markers <- read_excel(marker_dir,col_names = F)
  colnames(all_markers) <- c("cell_type","marker")
  markers.all <- c()
  for( i in 1:nrow(all_markers)){ ## i<-1
    cell_type <- all_markers[i,1]
    markers <- all_markers[i,2]
    markers <- unlist(strsplit(markers$marker, sep))
    markers <- gsub(" ","",markers)
    markers.all <- c(markers.all,markers) 
  }
  
  DotPlot(object = subset(obj,subset=cluV8=="1:Macro"), features = markers.all, assay = "SCT", group.by = "group", scale = TRUE, col.min = NULL, col.max = NULL)+
    scale_colour_gradient2(low = "#000000", mid = "orange", high = "red")+
    scale_size(range = c(1,10)) +
    labs(x=NULL,y=NULL,fill  = "Avg.exp",title = "Genes related to proinflammatory in macrophage")+
    scale_y_discrete(breaks=c("H","COE","CO"),labels=c("non-COVID","COVID-E","COVID-A"))+
    theme(axis.text.x = element_text(face = "bold",color="Black", size=12,hjust=1,vjust = 0.5,angle=90),
          axis.text.y = element_text(face = "bold",color="Black", size=12,hjust=1,vjust = 0.5,angle=0),
          legend.direction ="vertical",legend.position = "right")+
    guides(size=guide_legend(ncol=1,title = "Percent Expressed"))+
    guides(color = guide_colorbar(title = "Average Expression"))
  
  
}

# dot plot for 2: Fibroblast and 3: SMCs using gene list of “pro-fibrotic genes” and “anti-fibrotic genes” in excel.
if(T){
  ## read senescence associated genes
  marker_dir <- file.path(raw_data,"pro_fibrotic genes.xlsx")
  sep <- ","
  library(readxl)
  all_markers <- read_excel(marker_dir,col_names = F)
  colnames(all_markers) <- c("cell_type","marker")
  markers.all <- c()
  for( i in 1:nrow(all_markers)){ ## i<-1
    cell_type <- all_markers[i,1]
    markers <- all_markers[i,2]
    markers <- unlist(strsplit(markers$marker, sep))
    markers <- gsub(" ","",markers)
    markers.all <- c(markers.all,markers) 
  }
  
  ###2: Fibroblast and 3: SMCs
  DotPlot(object = subset(obj,subset = cluV8=="2:Fibro"), features = markers.all, assay = "SCT", group.by = "group", scale = TRUE, col.min = NULL, col.max = NULL)+
    scale_colour_gradient2(low = "#000000", mid = "orange", high = "red")+
    scale_size(range = c(1,10)) +
    labs(x=NULL,y=NULL,fill  = "Avg.exp",title = "Genes related to profibrosis in fibroblast")+
    scale_y_discrete(breaks=c("H","COE","CO"),labels=c("non-COVID","COVID-E","COVID-A"))+
    theme(axis.text.x = element_text(face = "bold",color="Black", size=12,hjust=1,vjust = 0.5,angle=90),
          axis.text.y = element_text(face = "bold",color="Black", size=12,hjust=1,vjust = 0.5,angle=0),
          legend.direction ="vertical",legend.position = "right")+
    guides(size=guide_legend(ncol=1,title = "Percent Expressed"))+
    guides(color = guide_colorbar(title = "Average Expression"))
  
  
}

## dot plot for 10: AT1, 5: AT2, 11: Basal, and 4: endothelial cells using gene list of 
## “senescence associated genes” and “SASP genes” in excel.
if(T){
  ## read senescence associated genes
  marker_dir <- file.path(raw_data,"senescence_related_genes.xlsx")
  sep <- ","
  library(readxl)
  all_markers <- read_excel(marker_dir,col_names = F)
  colnames(all_markers) <- c("cell_type","marker")
  markers.all <- c()
  for( i in 1:nrow(all_markers)){ ## i<-1
    cell_type <- all_markers[i,1]
    markers <- all_markers[i,2]
    markers <- unlist(strsplit(markers$marker, sep))
    markers <- gsub(" ","",markers)
    markers.all <- c(markers.all,markers) 
  }
  
  DotPlot(object = subset(obj,subset=cluV8=="5:AT2"), features = markers.all, assay = "SCT", group.by = "group", scale = TRUE, col.min = NULL, col.max = NULL)+
    scale_colour_gradient2(low = "#000000", mid = "orange", high = "red")+
    scale_size(range = c(1,10)) +
    labs(x=NULL,y=NULL,fill  = "Avg.exp",title = "Genes related to senescence in AT2")+
    scale_y_discrete(breaks=c("H","COE","CO"),labels=c("non-COVID","COVID-E","COVID-A"))+
    theme(axis.text.x = element_text(face = "bold",color="Black", size=12,hjust=1,vjust = 0.5,angle=90),
          axis.text.y = element_text(face = "bold",color="Black", size=12,hjust=1,vjust = 0.5,angle=0),
          legend.direction ="vertical",legend.position = "right")+
    guides(size=guide_legend(ncol=1,title = "Percent Expressed"))+
    guides(color = guide_colorbar(title = "Average Expression"))
  
  ## read SASP associated genes
  marker_dir <- file.path(raw_data,"SASP_related_genes.xlsx")
  sep <- ","
  library(readxl)
  all_markers <- read_excel(marker_dir,col_names = F)
  colnames(all_markers) <- c("cell_type","marker")
  markers.all <- c()
  for( i in 1:nrow(all_markers)){ ## i<-1
    cell_type <- all_markers[i,1]
    markers <- all_markers[i,2]
    markers <- unlist(strsplit(markers$marker, sep))
    markers <- gsub(" ","",markers)
    markers.all <- c(markers.all,markers) 
  }
  
  DotPlot(object = subset(obj,subset=cluV8=="5:AT2"), features = markers.all, assay = "SCT", group.by = "group", scale = TRUE, col.min = NULL, col.max = NULL)+
    scale_colour_gradient2(low = "#000000", mid = "orange", high = "red")+
    scale_size(range = c(1,10)) +
    labs(x=NULL,y=NULL,fill  = "Avg.exp",title = "Genes related to SASP in AT2")+
    scale_y_discrete(breaks=c("H","COE","CO"),labels=c("non-COVID","COVID-E","COVID-A"))+
    theme(axis.text.x = element_text(face = "bold",color="Black", size=12,hjust=1,vjust = 0.5,angle=90),
          axis.text.y = element_text(face = "bold",color="Black", size=12,hjust=1,vjust = 0.5,angle=0),
          legend.direction ="vertical",legend.position = "right")+
    guides(size=guide_legend(ncol=1,title = "Percent Expressed"))+
    guides(color = guide_colorbar(title = "Average Expression"))
  
  
  
}

### 8.	Generate images of each patient and each slide with all clusters.
### generate all images
if(T){
  ### functions used
  if(T){
    library(Seurat)
    library(data.table)
    library(sf)
    
    PolygonCropFOV <- function(fov, polygon, pyxel_size = 0.2125){
      #' Crop an FOV using a polygon
      #' 
      #' @description This function crops an FOV including centroids, segmentations and molecules.
      #'
      #' The FOV must contain centroids as those are used to determine if the cells are in the polygon or not. 
      #
      #' @param fov fov Object The FOV to be cropped.
      #' @param polygon sf POLYGON object.
      #' @param pixel_size ratio between the FOV coordinates and the sf Objects coordinates.
      #' @usage PolygonCropFOV(polygon, pyxel_size = 0.2125)
      #' @return A new cropped FOV object.
      #' @details If no cells are found the unchanged fov is returned with a Warning.
      #' @note TODO: if no centroids are present try to use segmentation or fail with an Error.
      
      # Use the centroids to identify the cells that should be kept
      message("Start processing centroids...")
      polygon <- st_set_crs(polygon, NA)
      centroids <- st_as_sfc(fov$centroids, forceMulti = FALSE)
      centroids <- st_set_crs(centroids, NA)
      centroids <- st_combine(centroids)
      centroids <- centroids / pyxel_size
      centroids <- st_cast(centroids, "POINT")
      contained <- st_contains_properly(polygon, centroids, sparse = FALSE)
      cells_to_keep <- Cells(fov)[contained]
      message("Done processing centroids...")
      
      # If no cells are found within the polygon, return the FOV unchanged with a warning
      if (!any(contained)){
        warning("No cells, found inside the polygon. Returnin g original object")
        return(fov)
      }
      
      # Build a new molecule object with only the molecules inside the polygon
      if (!is.null(fov$molecule)){
        message("Start processing molecules...")
        molecules <- rbindlist(lapply(fov$molecule, function(x){
          as.data.table(x@coords)
        }), use.names = T, idcol = "gene")
        mols <- st_multipoint(as.matrix(molecules[, .(x, y)]))
        mols <- mols / pyxel_size
        mols <- st_cast(st_sfc(mols), "POINT")
        mols <- st_set_crs(mols, NA)
        contained <- st_contains_properly(polygon, mols, sparse = FALSE)
        molecules <- molecules[as.logical(contained)]
        molecules <- CreateMolecules(as.data.frame(molecules))
        message("Done processing molecules...")
      }
      
      # Subset the input FOV to create a new FOV
      message("Preparing new fov...")
      new_fov <- subset(x = fov, cells = cells_to_keep) # This subsets the centroids and segmentation
      if (!is.null(fov$molecule)){
        new_fov[["molecules"]] <- molecules # Replace the molecules with the cropped ones
        new_fov[["molecule"]] <- molecules
      }
      message("Done cropping fov.")
      
      return(new_fov)
    }
    
    PolygonCropSeurat <- function(object, polygons, pyxel_size = 0.2125){
      #' Crop a Seurat Object using a list of Polygons
      #' 
      #' @description This function crops all FOVs included in the list, other FOVs included unchanged
      #' 
      #' The FOVs to be cropped must contain centroids as those are used to determine if the cells are in the polygon or not. 
      #
      #' @param object Seurat object to be cropped
      #' @param polygons a named list of sf POLYGON object, the names must match those of the FOVs to be cropped.
      #' @param pixel_size ratio between the FOV coordinates and the sf Objects coordinates.
      #' @usage PolygonCropSeurat(object polygons, pyxel_size = 0.2125)
      #' @return A new Seurat object, including the retained cells for the cropped FOVs and all cells for the other FOVs.
      #' @details It relies on `PolygonCropFOV` so if a cropped FOV would contain no cells, the uncropped FOV is included instead
      
      # Crop all FOVs for which a polygon has been provided
      message("Start cropping fovs...")
      cropped_fovs <- lapply(names(polygons), function(name){
        fov <- object@images[[name]]
        fov <- PolygonCropFOV(fov, polygons[[name]], pyxel_size)
        message("Done: ", name)
        return(fov)
      })
      names(cropped_fovs) <- names(polygons)
      message("Done cropping fovs...")
      
      # Make a list of all FOVs, cropped and non-cropped
      message("Gathering fovs...")
      untouched_fovs <- object@images[!(names(object@images) %in% names(cropped_fovs))]
      all_fovs <- c(untouched_fovs, cropped_fovs)
      message("Found fovs: ", names(all_fovs))
      
      # Get the cells to keep and make a subset of the object
      message("Making new object fovs...")
      cells_to_keep <- unlist(sapply(all_fovs, Cells))
      new_object <- subset(object, cells = cells_to_keep)
      new_object@images <- all_fovs # Replace the old FOVs
      message("Done cropping Seurat object")
      
      return(new_object)
    }
    
  }
  
  web.image <- file.path(raw_data,"image")
  if(!dir.exists(web.image)){dir.create(web.image,recursive = T)}
  
  all.genes <- row.names(obj)
  sysctl.pattern <- "^SystemControl"
  neg.pattern <- "^Negative"
  all.genes <- all.genes[which(! all.genes %in% c(grep(sysctl.pattern, all.genes, value=T), grep(neg.pattern, all.genes, value=T) ))]
  
  if(T){
    ### on slides, patients, fovs
    fov.image <- web.image
    for( IS in unique(obj@meta.data$slide)){ # IS <- "S7280.4"
      
      ### sub by slide ID for expression
      IS.sub <- subset_opt(obj, subset = slide == IS)
      # IS.sub <- subset(obj, subset = slide == IS)
      Idents(IS.sub) <- factor(IS.sub@meta.data$cluV8)
      
      col <- c("#66C5CC","#F89C74","#DCB0F2","#87C55F","#2F8AC4","#764E9F", "#DFFF00",
               "#0B775E","#9A8822","#7fcdbb","#dd1c77","#4D61C6","#0CD33F", "#696969" )
      
      ### sub by slide ID for segement
      options(future.globals.maxSize=80*1024^3)
      seg.xmin <- min(IS.sub$CenterX_global_px)
      seg.xmax <- max(IS.sub$CenterX_global_px)
      seg.ymin <- min(IS.sub$CenterY_global_px)
      seg.ymax <- max(IS.sub$CenterY_global_px)
      cropped.coords <- Crop(IS.sub[[unique(IS.sub$slide)]], x = c(seg.xmin, seg.xmax), y = c(seg.ymin, seg.ymax), coords = "tissue")
      IS.sub[["zoom1"]] <- cropped.coords
      DefaultBoundary(IS.sub[["zoom1"]]) <- "segmentation"
      
      ### sub extract the fov for patient
      for(IP in unique(IS.sub@meta.data$patient)){ # IP <- "H06"
        ### create the folder for patient
        fov.image.slide.patient <- file.path(fov.image,IP)
        if(!dir.exists(fov.image.slide.patient)){dir.create(fov.image.slide.patient,recursive = T)}
        
        ## subset
        IP.sub <- subset_opt(IS.sub, subset = patient == IP)
        Idents(IP.sub) <- factor(IP.sub@meta.data$cluV8)
        
        ### sub by slide ID for segement
        seg.xmin <- min(IP.sub$CenterX_global_px)
        seg.xmax <- max(IP.sub$CenterX_global_px)
        seg.ymin <- min(IP.sub$CenterY_global_px)
        seg.ymax <- max(IP.sub$CenterY_global_px)
        cropped.coords <- Crop(IP.sub[[unique(IP.sub$slide)]], x = c(seg.xmin, seg.xmax), y = c(seg.ymin, seg.ymax), coords = "tissue")
        IP.sub[["zoom1"]] <- cropped.coords
        DefaultBoundary(IP.sub[["zoom1"]]) <- "segmentation"
        
        ### sub exract the fov for each fov
        for(IF in unique(IP.sub@meta.data$fov)){ ## IF <- 36
          ### subset
          IF.sub <- subset_opt(IP.sub, subset = fov == IF)
          Idents(IF.sub) <- factor(IF.sub@meta.data$cluV8)
          
          ### fake to zoom in the fov
          seg.xmin <- min(IF.sub$CenterX_global_px)
          seg.xmax <- max(IF.sub$CenterX_global_px)
          seg.ymin <- min(IF.sub$CenterY_global_px)
          seg.ymax <- max(IF.sub$CenterY_global_px)
          cropped.coords <- Crop(IF.sub[[unique(IF.sub$slide)]], x = c(seg.xmin, seg.xmax), y = c(seg.ymin, seg.ymax), coords = "tissue")
          IF.sub[["zoom1"]] <- cropped.coords
          DefaultBoundary(IF.sub[["zoom1"]]) <- "segmentation"
          
          g <- ImageDimPlot(IF.sub, fov = "zoom1", cols = col, alpha = 0.6, crop=T, axes=T, dark.background=F,
                            mols.size = 1.5, nmols = 20000, border.color = NA, coord.fixed = T,size=1,mols.alpha=1)+
            theme_bw()+ theme(panel.grid.minor = element_blank(),  panel.grid.major = element_blank() )+
            scale_fill_manual(breaks = c("0:LC-1","1:Macro","2:Fibro","3:SMC","4:Endo","5:AT2","6:Mono",      
                                         "7:Plasma","8:RBC","9:LC-2","10:AT1","11:Basal","12:Nutrophil","13:Unknown"),
                              values = c("#66C5CC","#F89C74","#DCB0F2","#87C55F","#2F8AC4","#764E9F", "#DFFF00",
                                         "#0B775E","#9A8822","#7fcdbb","#dd1c77","#4D61C6","#0CD33F", "#696969" ))
          print(g)
          ggsave(file.path(fov.image.slide.patient,paste0(unique(IF.sub$patient),".",IF,".png")),width = 8,height = 7,dpi = 300,plot = g)
          # 
          
          for(i in all.genes){ # i <- "CXCR6"
              g <- ImageDimPlot(IF.sub, fov = "zoom1", cols = col, alpha = 0.3, molecules = i,crop=T, axes=T, dark.background=F,
                                mols.cols = "red",mols.size = 1.5, nmols = 20000, border.color = NA, coord.fixed = T,size=1,mols.alpha=1)+
                  theme_bw()+ theme(panel.grid.minor = element_blank(),  panel.grid.major = element_blank() )+
                  guides(fill = guide_legend(override.aes = list(alpha = 0.3)))+
                  scale_fill_manual(breaks = c("0:LC-1","1:Macro","2:Fibro","3:SMC","4:Endo","5:AT2","6:Mono",
                                               "7:Plasma","8:RBC","9:LC-2","10:AT1","11:Basal","12:Nutrophil","13:Unknown"),
                                    values = c("#66C5CC","#F89C74","#DCB0F2","#87C55F","#2F8AC4","#764E9F", "#DFFF00",
                                               "#0B775E","#9A8822","#7fcdbb","#dd1c77","#4D61C6","#0CD33F", "#696969" ))
              print(g)
              ggsave(file.path(fov.image.slide.patient,paste0(unique(IF.sub$patient),".",IF,".",gsub("\\/",".",i),".png")),width = 8,height = 7,dpi = 300,plot = g)

          }
        }
      }
    }
  }
}

### Cell number of each FOV
if(T){
  library(dplyr)
  library(tidyr)
  df <- obj@meta.data %>% as.data.frame() %>% subset(., select = c("group","patient","fov","cluV8")) %>% 
    group_by(across(all_of(c("group","patient","fov","cluV8")))) %>% dplyr::summarise(cell_count = n()) %>% 
    tidyr::spread(key = "cluV8", value = "cell_count") %>% as.data.frame()
  df[is.na(df)] <- 0
  write_xlsx(df,file.path(raw_data,"06.cell number for each fov.xlsx"),col_names = T)
  
}

### for the number of macrophage, monocyte, Lymphocyte, and Nutrophil
if(T){
  ## read the cell number file
  cell.Num.dir <- file.path(raw_data,"06.cell number for each fov.xlsx")
  sep <- ","
  library(readxl)
  cell.Num <- read_excel(cell.Num.dir,col_names = T)
  
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  library(ggsignif)
  library(ggbeeswarm)
  
  ### for Macrophage
  cell.Num.macro <- cell.Num[,c("group","1:Macro")]
  cell.Num.macro <- tidyr::gather(cell.Num.macro,key="1:Macro", value="value", -"group")
  cell.Num.macro$Group <- factor(cell.Num.macro$group,levels = c("H","COE","CO"))
  
  ##calculate the p-value
  group.pvalue.CO.vs.H <- t.test(cell.Num.macro$value[cell.Num.macro$Group == "H"], 
                                 cell.Num.macro$value[cell.Num.macro$Group == "CO"])
  
  group.pvalue.COE.vs.CO <- t.test(cell.Num.macro$value[cell.Num.macro$Group == "CO"], 
                                   cell.Num.macro$value[cell.Num.macro$Group == "COE"])
  
  group.pvalue.COE.vs.H <- t.test(cell.Num.macro$value[cell.Num.macro$Group == "H"], 
                                  cell.Num.macro$value[cell.Num.macro$Group == "COE"])
  
  g1 <-  ggplot(cell.Num.macro, aes(x=Group, y=value, group=Group, color=Group)) +
    stat_summary(fun = "mean", geom = "crossbar",size = 0.2,width=0.4,color = "black")+
    geom_beeswarm(cex=0.9,dodge.width=0.1,size=1,alpha =.6) +
    geom_violin(alpha = 0.5)+
    geom_signif(annotations = c("p-value = 9.83e-05","p-value = 8.02e-04","p-value = 0.23"), y_position = c(800,700,600), xmin = c(1,2,1),
                xmax = c(3,3,2), tip_length = c(0.02,0.02,0.02,0.02,0.02,0.02),color="black",size=0.5) +
    labs(x="", y="# of Macrophage",title = "# of Macropahge")+
    scale_x_discrete(breaks=c("H","COE","CO"), labels=c("non-COVID", "COVID-E", "COVID-A"))+
    scale_color_manual(values=c("#b8edc6", "#82b2ff", "#ffada7"))+
    theme_bw()+
    theme(axis.text.x = element_text(face = "bold",color="Black", size=12,hjust=0.5,vjust = 0.5,angle=0),
          axis.text.y = element_text(face = "bold",color="Black", size=8,hjust=1,vjust = 0.5,angle=0),
          axis.title.x = element_text(face="bold", color ="Black", size = 10, angle = 0, hjust = .5, vjust = 0.5),
          axis.title.y = element_text(face="bold",color ="Black", size = 12, angle = 90, hjust = 0.5, vjust = 1.5),
          legend.position="none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          strip.text = element_text(face="bold", color ="Black", size = 12)
    )
  
  ### for Monocyte
  cell.Num.mono <- cell.Num[,c("group","6:Mono")]
  cell.Num.mono <- tidyr::gather(cell.Num.mono,key="6:Mono", value="value", -"group")
  cell.Num.mono$Group <- factor(cell.Num.mono$group,levels = c("H","COE","CO"))
  
  ##calculate the p-value
  group.pvalue.CO.vs.H <- t.test(cell.Num.mono$value[cell.Num.mono$Group == "H"], 
                                 cell.Num.mono$value[cell.Num.mono$Group == "CO"])
  
  group.pvalue.COE.vs.CO <- t.test(cell.Num.mono$value[cell.Num.mono$Group == "CO"], 
                                   cell.Num.mono$value[cell.Num.mono$Group == "COE"])
  
  group.pvalue.COE.vs.H <- t.test(cell.Num.mono$value[cell.Num.mono$Group == "H"], 
                                  cell.Num.mono$value[cell.Num.mono$Group == "COE"])
  
  g2 <-  ggplot(cell.Num.mono, aes(x=Group, y=value, group=Group, color=Group)) +
    stat_summary(fun = "mean", geom = "crossbar",size = 0.2,width=0.4,color = "black")+
    geom_beeswarm(cex=0.9,dodge.width=0.1,size=1,alpha =.6) +
    geom_violin(alpha = 0.5)+
    geom_signif(annotations = c("p-value = 1.12e-06","p-value = 0.03","p-value = 0.12"), y_position = c(500,380,420), xmin = c(1,1,2),
                xmax = c(3,2,3), tip_length = c(0.02,0.02,0.02,0.02,0.02,0.02),color="black",size=0.5) +
    labs(x="", y="# of Monocyte",title = "# of Monocyte")+
    scale_x_discrete(breaks=c("H","COE","CO"), labels=c("non-COVID", "COVID-E", "COVID-A"))+
    scale_color_manual(values=c("#b8edc6", "#82b2ff", "#ffada7"))+
    theme_bw()+
    theme(axis.text.x = element_text(face = "bold",color="Black", size=12,hjust=0.5,vjust = 0.5,angle=0),
          axis.text.y = element_text(face = "bold",color="Black", size=8,hjust=1,vjust = 0.5,angle=0),
          axis.title.x = element_text(face="bold", color ="Black", size = 10, angle = 0, hjust = .5, vjust = 0.5),
          axis.title.y = element_text(face="bold",color ="Black", size = 12, angle = 90, hjust = 0.5, vjust = 1.5),
          legend.position="none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          strip.text = element_text(face="bold", color ="Black", size = 12)
    )
  
  
  ### for neutriphils
  cell.Num.neutriphils <- cell.Num[,c("group","12:Nutrophil")]
  cell.Num.neutriphils <- tidyr::gather(cell.Num.neutriphils,key="12:Nutrophil", value="value", -"group")
  cell.Num.neutriphils$Group <- factor(cell.Num.neutriphils$group,levels = c("H","COE","CO"))
  
  ##calculate the p-value
  group.pvalue.CO.vs.H <- t.test(cell.Num.neutriphils$value[cell.Num.neutriphils$Group == "H"], 
                                 cell.Num.neutriphils$value[cell.Num.neutriphils$Group == "CO"])
  
  group.pvalue.COE.vs.CO <- t.test(cell.Num.neutriphils$value[cell.Num.neutriphils$Group == "CO"], 
                                   cell.Num.neutriphils$value[cell.Num.neutriphils$Group == "COE"])
  
  group.pvalue.COE.vs.H <- t.test(cell.Num.neutriphils$value[cell.Num.neutriphils$Group == "H"], 
                                  cell.Num.neutriphils$value[cell.Num.neutriphils$Group == "COE"])
  
  g3 <-  ggplot(cell.Num.neutriphils, aes(x=Group, y=value, group=Group, color=Group)) +
    stat_summary(fun = "mean", geom = "crossbar",size = 0.2,width=0.4,color = "black")+
    geom_beeswarm(cex=0.5,dodge.width=0.1,size=1,alpha =.6) +
    geom_violin(alpha = 0.5)+
    geom_signif(annotations = c("p-value = 6.92e-05","p-value = 2.71e-11","p-value = 0.23"), y_position = c(600,500,400), xmin = c(1,2,1),
                xmax = c(3,3,2), tip_length = c(0.02,0.02,0.02,0.02,0.02,0.02),color="black",size=0.5) +
    labs(x="", y="# of Neutrophil",title = "# of Neutrophil")+
    scale_x_discrete(breaks=c("H","COE","CO"), labels=c("non-COVID", "COVID-E", "COVID-A"))+
    scale_color_manual(values=c("#b8edc6", "#82b2ff", "#ffada7"))+
    coord_flip()+
    theme_bw()+
    theme(axis.text.x = element_text(face = "bold",color="Black", size=12,hjust=0.5,vjust = 0.5,angle=0),
          axis.text.y = element_text(face = "bold",color="Black", size=8,hjust=1,vjust = 0.5,angle=0),
          axis.title.x = element_text(face="bold", color ="Black", size = 10, angle = 0, hjust = .5, vjust = 0.5),
          axis.title.y = element_text(face="bold",color ="Black", size = 12, angle = 90, hjust = 0.5, vjust = 1.5),
          legend.position="none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          strip.text = element_text(face="bold", color ="Black", size = 12)
    )
  
  ### for LC
  cell.Num.LC <- cell.Num[,c("group","0:LC-1","9:LC-2")]
  cell.Num.LC$LC <- apply(cell.Num.LC[,c("0:LC-1","9:LC-2")],1,sum)
  cell.Num.LC <- cell.Num.LC[,c("group","LC")]
  cell.Num.LC <- tidyr::gather(cell.Num.LC,key="LC", value="value", -"group")
  cell.Num.LC$Group <- factor(cell.Num.LC$group,levels = c("H","COE","CO"))
  
  ##calculate the p-value
  group.pvalue.CO.vs.H <- t.test(cell.Num.LC$value[cell.Num.LC$Group == "H"], 
                                 cell.Num.LC$value[cell.Num.LC$Group == "CO"])
  
  group.pvalue.COE.vs.CO <- t.test(cell.Num.LC$value[cell.Num.LC$Group == "CO"], 
                                   cell.Num.LC$value[cell.Num.LC$Group == "COE"])
  
  group.pvalue.COE.vs.H <- t.test(cell.Num.LC$value[cell.Num.LC$Group == "H"], 
                                  cell.Num.LC$value[cell.Num.LC$Group == "COE"])
  
  g4 <-  ggplot(cell.Num.LC, aes(x=Group, y=value, group=Group, color=Group)) +
    stat_summary(fun = "mean", geom = "crossbar",size = 0.2,width=0.4,color = "black")+
    geom_beeswarm(cex=0.5,dodge.width=0.1,size=1,alpha =.6) +
    geom_violin(alpha = 0.5)+
    geom_signif(annotations = c("p-value = 0.70","p-value = 0.73","p-value = 0.50"), y_position = c(900,1100,1300), xmin = c(1,2,1),
                xmax = c(2,3,3), tip_length = c(0.02,0.02,0.02,0.02,0.02,0.02),color="black",size=0.5) +
    labs(x="", y="# of Lymphocyte",title = "# of Lymphocyte")+
    scale_x_discrete(breaks=c("H","COE","CO"), labels=c("non-COVID", "COVID-E", "COVID-A"))+
    scale_color_manual(values=c("#b8edc6", "#82b2ff", "#ffada7"))+
    coord_flip()+
    theme_bw()+
    theme(axis.text.x = element_text(face = "bold",color="Black", size=12,hjust=0.5,vjust = 0.5,angle=0),
          axis.text.y = element_text(face = "bold",color="Black", size=8,hjust=1,vjust = 0.5,angle=0),
          axis.title.x = element_text(face="bold", color ="Black", size = 10, angle = 0, hjust = .5, vjust = 0.5),
          axis.title.y = element_text(face="bold",color ="Black", size = 12, angle = 90, hjust = 0.5, vjust = 1.5),
          legend.position="none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          strip.text = element_text(face="bold", color ="Black", size = 12)
    )
  
  pdf(file.path(raw_data,"11.cell number comparison1.pdf"),width = 12,height = 6)
  print(g1|g2 )
  dev.off()
  
  pdf(file.path(raw_data,"11.cell number comparison2.pdf"),width = 12,height = 5)
  print(g3|g4 )
  dev.off()
  
  
}

# add the quantification for IGFBP3/6 in AT2
if(T){
  # IGFBP6
  ### FOV:
  # H02(39), H03(23), H05(5), H06(21), H08(38), H09(55)
  # CO01(20), CO02(86), CO03(8), CO05(59), CP07(69), CO12(94), CO16(74), CO19(3), CO20(16), CO21(55)
  # COE01(107), COE02(87), COE04(43), COE05(91), COE06(27), COE07(49), 
  # calculate the cell number expressed IGFBP6 in each fov and compare among groups
  obj.igfbp6 <- subset_opt(obj, subset=  (patient == "CO01" & fov %in% c(20)) | (patient == "CO02" & fov %in% c(86)) | 
                             (patient == "CO03" & fov %in% c(8)) | (patient == "CO05" & fov %in% c(59)) |
                             (patient == "CP07" & fov %in% c(69)) | (patient == "CO12" & fov %in% c(94)) | 
                             (patient == "CO16" & fov %in% c(74)) | (patient == "CO19" & fov %in% c(3)) | 
                             (patient == "CO20" & fov %in% c(16)) | (patient == "CO21" & fov %in% c(55)) |
                             (patient == "COE01" & fov %in% c(107)) | (patient == "COE02" & fov %in% c(87)) |
                             (patient == "COE04" & fov %in% c(43)) | (patient == "COE05" & fov %in% c(91)) |
                             (patient == "COE06" & fov %in% c(27)) | (patient == "COE07" & fov %in% c(49)) |
                             (patient == "H02" & fov %in% c(39)) | (patient == "H03" & fov %in% c(23)) |
                             (patient == "H05" & fov %in% c(5)) | (patient == "H06" & fov %in% c(21)) |
                             (patient == "H08" & fov %in% c(38)) | (patient == "H09" & fov %in% c(55))  )
  
  obj.igfbp6$igfbp6.exp <- ifelse(obj.igfbp6@assays$SCT$data["IGFBP6", ] >= 0, 1, 0)
  
  obj.igfbp6$group.sample.fov <- paste0(obj.igfbp6$group,"_",obj.igfbp6$patient,"_",obj.igfbp6$fov)
  
  igfbp6_counts <- obj.igfbp6@meta.data %>%
    group_by(cell_type = obj.igfbp6$cluV8, sample = obj.igfbp6$group.sample.fov) %>%
    summarise(igfbp6_positive_cells = sum(igfbp6.exp),
              total_cells = n(),
              proportion_igfbp6_positive = igfbp6_positive_cells / total_cells) %>%
    ungroup()
  
  igfbp6_counts.at2 <- igfbp6_counts[which(igfbp6_counts$cell_type == "5:AT2"),]
  
  library(dplyr)
  library(tidyr)
  igfbp6_counts.at2 <- igfbp6_counts.at2 %>% separate(sample, c('group', "patient","fov"),sep="_")
  
  
  #calculate the number of positive IGFBP6 AT2 cells
  ### for AT2
  cell.Num.AT2 <- igfbp6_counts.at2[,c("group","igfbp6_positive_cells")]
  # cell.Num.AT2 <- igfbp6_counts.at2[,c(2,4)]
  colnames(cell.Num.AT2) <- c("Group","value")
  cell.Num.AT2$Group <- factor(cell.Num.AT2$Group,levels = c("H","COE","CO"))
  cell.Num.AT2$value <-as.numeric(cell.Num.AT2$value)
  
  ##calculate the p-value
  group.pvalue.CO.vs.H <- t.test(cell.Num.AT2$value[cell.Num.AT2$Group == "H"], 
                                 cell.Num.AT2$value[cell.Num.AT2$Group == "CO"])
  
  group.pvalue.COE.vs.CO <- t.test(cell.Num.AT2$value[cell.Num.AT2$Group == "CO"], 
                                   cell.Num.AT2$value[cell.Num.AT2$Group == "COE"])
  
  group.pvalue.COE.vs.H <- t.test(cell.Num.AT2$value[cell.Num.AT2$Group == "H"], 
                                  cell.Num.AT2$value[cell.Num.AT2$Group == "COE"])
  library(ggbeeswarm)
  g2 <-  ggplot(cell.Num.AT2, aes(x=Group, y=value, group=Group, color=Group)) +
    # stat_summary(fun.data=data_summary,geom="errorbar",width=0.2,size=0.5,color = "black") +
    stat_summary(fun = "mean", geom = "crossbar",size = 0.2,width=0.4,color = "black")+
    # geom_dotplot(binaxis='y', stackdir='center',binwidth=0.06,stackratio=1.2, dotsize=0.8) +
    geom_beeswarm(cex=0.9,dodge.width=0.3,size=1,alpha =1) +
    geom_violin(alpha = 0.5)+
    geom_signif(annotations = c(paste0("p-value = ",format(round(group.pvalue.CO.vs.H$p.value,digits = 2),nsmall = 2)),
                                paste0("p-value = ",format(round(group.pvalue.COE.vs.CO$p.value,digits = 2),nsmall = 2)),
                                paste0("p-value = ",format(round(group.pvalue.COE.vs.H$p.value,digits = 2),nsmall = 2)) ), 
                y_position = c(260,290,320), 
                xmin = c(1,2,1),
                xmax = c(2,3,3), tip_length = c(0.0,0.0),color="black",size=0.5) +
    # coord_cartesian(ylim=c(9.5,15.5))+
    # scale_y_continuous(breaks = c(10,11,12,13,14,15)) +
    labs(x="", y="# of IGFBP6+ AT2 cells")+
    scale_x_discrete(breaks=c("H","COE","CO"), labels=c("non-COVID", "COVID-E", "COVID-A"))+
    scale_color_manual(values=c("#abd9b7", "#6996df", "#e17a74"))+
    theme_bw()+
    theme(axis.text.x = element_text(face = "bold",color="Black", size=12,hjust=0.5,vjust = 0.5,angle=0),
          axis.text.y = element_text(face = "bold",color="Black", size=8,hjust=1,vjust = 0.5,angle=0),
          axis.title.x = element_text(face="bold", color ="Black", size = 10, angle = 0, hjust = .5, vjust = 0.5),
          axis.title.y = element_text(face="bold",color ="Black", size = 12, angle = 90, hjust = 0.5, vjust = 1.5),
          # legend.position=c(0.88,0.75), legend.key = element_blank(), legend.background = element_blank(),
          # legend.key.size = unit(6, "pt"),
          # legend.text = element_text(color ="Black", size = 9, angle = 0, hjust = 0, vjust = 0.5),
          # legend.title = element_text(face="bold",color ="Black", size = 10, angle = 0, hjust = .5, vjust = .5),
          legend.position="none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          # axis.ticks=element_blank(),
          panel.border = element_blank(),
          strip.text = element_text(face="bold", color ="Black", size = 12)
    )
  print(g2)
  
  pdf(file.path(raw_data,"12.quantification of IGFBP6+ cell number.pdf"),width = 4,height = 6)
  print(g2)
  dev.off()
  
  
}

