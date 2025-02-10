#Count arm-level aneuploidy using set chromosome thresholds
library(data.table)
library(stringr)

log2_to_CNA_needed <- TRUE 

#note: 2019 CCLE data is more conservative in CNA calls, whereas 2012 data matches closer to TCGA reporting
percent_of_arm_threshold <- 66.7
input_CNAs <- "data_CCLEcna_2012_Nature_hg19.seg"
output_arm_CNAs <- "CCLE_arm_signed_aneuploidy.tsv"

#find the first and last gene locations
arm_edge_genes <- as.data.frame(fread("arm_edges.tsv"))

#keep only a threshold of SCNA, of 0.3
if(str_detect(input_CNAs,".seg")){
  input_CNAs <- as.data.frame(fread(input_CNAs))
  colnames(input_CNAs) <- c("sampleID","chr","start","end","markers","value")
  input_CNAs <- input_CNAs[which(input_CNAs$markers >= 100),]
  input_CNAs$markers <- NULL
  }else{
  input_CNAs <- as.data.frame(fread(input_CNAs))
}

input_CNAs <- input_CNAs[which( abs(input_CNAs$value) >= 0.3),]

if(log2_to_CNA_needed == TRUE){
  CNA_cutoff_range <- c(-1, -0.3, 0.3, 1)
  isLoss <- function(x){if(is.na(x)==FALSE & is.null(x)==FALSE){if(x<CNA_cutoff_range[2] & x>CNA_cutoff_range[1]){TRUE}else{FALSE}}else{NA}}
  isGain <- function(x){if(is.na(x)==FALSE & is.null(x)==FALSE){if(x<CNA_cutoff_range[4] & x>CNA_cutoff_range[3]){TRUE}else{FALSE}}else{NA}}
  isNone <- function(x){if(is.na(x)==FALSE & is.null(x)==FALSE){if(x<CNA_cutoff_range[3] & x>CNA_cutoff_range[2]){TRUE}else{FALSE}}else{NA}}
  isDeletion <- function(x){if(is.na(x)==FALSE & is.null(x)==FALSE){if(x<CNA_cutoff_range[1]){TRUE}else{FALSE}}else{NA}}
  isAmplification <- function(x){if(is.na(x)==FALSE & is.null(x)==FALSE){if(x>CNA_cutoff_range[4]){TRUE}else{FALSE}}else{NA}}
  
  input_CNAs$value <- input_CNAs$value*0 +
    as.integer(sapply(input_CNAs$value, isDeletion)*-2) +
    as.integer(sapply(input_CNAs$value, isLoss)*-1) +
    as.integer(sapply(input_CNAs$value, isGain)*1) +
    as.integer(sapply(input_CNAs$value, isAmplification)*2)
  
}

if(str_detect(input_CNAs$chr[1],"chr") == FALSE){
  input_CNAs$chr <- paste0("chr",input_CNAs$chr)
}

arm_edge_genes$arm_length <- arm_edge_genes$arm_end - arm_edge_genes$arm_start
input_CNAs <- input_CNAs[!input_CNAs$chr %in% c("chrX","chrY"),] #remove sex chromosome data
input_CNAs$chrp_length <- 0
input_CNAs$chrq_length <- 0

#convert SCNAs to arm level
chrom_label <- NULL
input_CNAs_chrom <- input_CNAs[1,]
input_CNAs_chrom_p <- input_CNAs[1,]
input_CNAs_chrom_q <- input_CNAs[1,]
input_CNAs_chrommixed <- input_CNAs[1,]


input_CNAs_arm_split <- lapply(unique(input_CNAs$chr), function(chrom_label){ #chrom_label <- "chr21"
  if(paste0(chrom_label,"p") %in% arm_edge_genes$chrom_arm){
  input_CNAs_chrom <- input_CNAs[which(input_CNAs$chr == chrom_label),]
  input_CNAs_chrom$chrp_only <- as.integer(input_CNAs_chrom$end <= arm_edge_genes$arm_end[which(arm_edge_genes$chrom_arm == paste0(chrom_label,"p"))])
  input_CNAs_chrom$chrq_only <- as.integer(input_CNAs_chrom$start >= arm_edge_genes$arm_start[which(arm_edge_genes$chrom_arm == paste0(chrom_label,"q"))])
  input_CNAs_chrom_p$chrp_length <- 0
  input_CNAs_chrom_p$chrq_length <- 0
  
  input_CNAs_chrommixed <- input_CNAs_chrom[which((input_CNAs_chrom$chrp_only + input_CNAs_chrom$chrq_only)==0),]
  input_CNAs_chrommixed$chrp_diff <- as.integer(input_CNAs_chrommixed$start - arm_edge_genes$arm_end[which(arm_edge_genes$chrom_arm == paste0(chrom_label,"p"))])
  input_CNAs_chrommixed$chrq_diff <- as.integer(input_CNAs_chrommixed$end - arm_edge_genes$arm_start[which(arm_edge_genes$chrom_arm == paste0(chrom_label,"q"))])
  input_CNAs_chrommixed$chrp_length <- -input_CNAs_chrommixed$chrp_diff
  input_CNAs_chrommixed$chrp_length[which(input_CNAs_chrommixed$chrp_diff > 0)] <- 0
  input_CNAs_chrommixed$chrq_length <- input_CNAs_chrommixed$chrq_diff
  input_CNAs_chrommixed$chrq_length[which(input_CNAs_chrommixed$chrq_diff < 0)] <- 0
  input_CNAs_chrommixed <- input_CNAs_chrommixed[,colnames(input_CNAs)]
  
  input_CNAs_chrom_p <- input_CNAs_chrom[which(input_CNAs_chrom$chrp_only == 1),]
  input_CNAs_chrom_p$chrp_length <- input_CNAs_chrom_p$end - input_CNAs_chrom_p$start
  input_CNAs_chrom_p <- rbind(input_CNAs_chrom_p[,colnames(input_CNAs)], input_CNAs_chrommixed[,colnames(input_CNAs)])
  input_CNAs_chrom_p <- input_CNAs_chrom_p[which(input_CNAs_chrom_p$chrp_length > 0),]
  input_CNAs_chrom_p <- input_CNAs_chrom_p[,colnames(input_CNAs)]
  
  input_CNAs_chrom_q <- input_CNAs_chrom[which(input_CNAs_chrom$chrq_only == 1),]
  input_CNAs_chrom_q$chrq_length <- input_CNAs_chrom_q$end - input_CNAs_chrom_q$start
  input_CNAs_chrom_q <- rbind(input_CNAs_chrom_q[,colnames(input_CNAs)], input_CNAs_chrommixed[,colnames(input_CNAs)])
  input_CNAs_chrom_q <- input_CNAs_chrom_q[which(input_CNAs_chrom_q$chrq_length > 0),]
  input_CNAs_chrom_q <- input_CNAs_chrom_q[,colnames(input_CNAs)]
  
  input_CNAs_chrom <- rbind(input_CNAs_chrom_p,input_CNAs_chrom_q, input_CNAs_chrommixed)
  input_CNAs_chrom <- input_CNAs_chrom[,colnames(input_CNAs)]
  input_CNAs_chrom <- input_CNAs_chrom[!duplicated(input_CNAs_chrom),]
  
  } else { #if no p arm
    input_CNAs_chrom <- input_CNAs[which(input_CNAs$chr == chrom_label),]
    input_CNAs_chrom$chrp_length <- 0
    input_CNAs_chrom$start[which(input_CNAs_chrom$start <  arm_edge_genes$arm_start[which(arm_edge_genes$chrom_arm == paste0(chrom_label,"q"))])] <- as.integer(arm_edge_genes$arm_start[which(arm_edge_genes$chrom_arm == paste0(chrom_label,"q"))])
    input_CNAs_chrom$chrq_length <- input_CNAs_chrom$end - input_CNAs_chrom$start
    input_CNAs_chrom <- input_CNAs_chrom[which(input_CNAs_chrom$chrq_length > 0),] #filters regions w/o genes or in p/cen region
  }
  
  
  return(input_CNAs_chrom)
})


input_CNAs <- do.call("rbind", input_CNAs_arm_split)
input_CNAs <- input_CNAs[which(input_CNAs$chrp_length + input_CNAs$chrq_length != 0),]


#####All alterations#####
output_matrix <- matrix(0, nrow = length(unique(input_CNAs$sampleID)), ncol = nrow(arm_edge_genes))
row.names(output_matrix) <- unique(input_CNAs$sampleID)
colnames(output_matrix) <- arm_edge_genes$chrom_arm
sample_names <- unique(input_CNAs$sampleID)

input_CNAs_sample <- input_CNAs[1,]
for(sample_name in sample_names){ #sample_name <- sample_names[3]
  input_CNAs_sample <-  input_CNAs[which( input_CNAs$sampleID == sample_name),]
  for(chrom_num in colnames(output_matrix)){ #chrom_num <- "chr5p"
    if(str_detect(chrom_num,"p") == TRUE){
    output_matrix[sample_name, chrom_num] <- sum(input_CNAs_sample$chrp_length[
      which(paste0(input_CNAs_sample$chr,"p") == chrom_num)])
    } else { 
      if(str_detect(chrom_num,"q") == TRUE){ #QC to ensure q arm exists
        output_matrix[sample_name, chrom_num] <- sum(input_CNAs_sample$chrq_length[
          which(paste0(input_CNAs_sample$chr,"q") == chrom_num)])
      }
      }
    }
}

length_divisors <- arm_edge_genes$arm_end -   arm_edge_genes$arm_start
names(length_divisors) <- arm_edge_genes$chrom_arm
for(chrom_arm in colnames(output_matrix)){
  output_matrix[,chrom_arm] <- round(output_matrix[,chrom_arm] / length_divisors[chrom_arm] * 100 , digits = 3)
}

#some samples overlap regions, ensure none are > 100%
output_matrix[which(output_matrix > 100)] <- 100



#####Loss alterations#####

output_matrix_loss <- matrix(0, nrow = length(unique(input_CNAs$sampleID)), ncol = nrow(arm_edge_genes))
row.names(output_matrix_loss) <- unique(input_CNAs$sampleID)
colnames(output_matrix_loss) <- arm_edge_genes$chrom_arm
sample_names <- unique(input_CNAs$sampleID)

input_CNAs_sample <- input_CNAs[1,]
for(sample_name in sample_names){ #sample_name <- sample_names[3]
  input_CNAs_sample <-  input_CNAs[which( input_CNAs$sampleID == sample_name),]
  input_CNAs_sample <-  input_CNAs_sample[which(input_CNAs_sample$value < 0 ),]
  input_CNAs_sample <-  input_CNAs_sample[!is.na(input_CNAs_sample$start),]
  for(chrom_num in colnames(output_matrix_loss)){ #chrom_num <- "chr5p"
    if(str_detect(chrom_num,"p") == TRUE){
      output_matrix_loss[sample_name, chrom_num] <- sum(input_CNAs_sample$chrp_length[
        which(paste0(input_CNAs_sample$chr,"p") == chrom_num)])
    } else { 
      if(str_detect(chrom_num,"q") == TRUE){ #QC to ensure q arm exists
        output_matrix_loss[sample_name, chrom_num] <- sum(input_CNAs_sample$chrq_length[
          which(paste0(input_CNAs_sample$chr,"q") == chrom_num)])
      }
    }
  }
}

for(chrom_arm in colnames(output_matrix_loss)){
  output_matrix_loss[,chrom_arm] <- round(output_matrix_loss[,chrom_arm] / length_divisors[chrom_arm] * 100 , digits = 3)
}

#if samples overlap regions, ensure none are > 100%
output_matrix_loss[which(output_matrix_loss > 100)] <- 100



#####Gain alterations#####
output_matrix_gain <- matrix(0, nrow = length(unique(input_CNAs$sampleID)), ncol = nrow(arm_edge_genes))
row.names(output_matrix_gain) <- unique(input_CNAs$sampleID)
colnames(output_matrix_gain) <- arm_edge_genes$chrom_arm
sample_names <- unique(input_CNAs$sampleID)

input_CNAs_sample <- input_CNAs[1,]
for(sample_name in sample_names){ #sample_name <- sample_names[1]
  input_CNAs_sample <-  input_CNAs[which( input_CNAs$sampleID == sample_name),]
  input_CNAs_sample <-  input_CNAs_sample[which(input_CNAs_sample$value > 0 ),]
  input_CNAs_sample <-  input_CNAs_sample[!is.na(input_CNAs_sample$start),]
  for(chrom_num in colnames(output_matrix_gain)){ #chrom_num <- "chr5p"
    if(str_detect(chrom_num,"p") == TRUE){
      output_matrix_gain[sample_name, chrom_num] <- sum(input_CNAs_sample$chrp_length[
        which(paste0(input_CNAs_sample$chr,"p") == chrom_num)])
    } else { 
      if(str_detect(chrom_num,"q") == TRUE){ #QC to ensure q arm exists
        output_matrix_gain[sample_name, chrom_num] <- sum(input_CNAs_sample$chrq_length[
          which(paste0(input_CNAs_sample$chr,"q") == chrom_num)])
      }
    }
  }
}

for(chrom_arm in colnames(output_matrix_gain)){
  output_matrix_gain[,chrom_arm] <- round(output_matrix_gain[,chrom_arm] / length_divisors[chrom_arm] * 100 , digits = 3)
}

#some samples overlap regions, ensure none are > 100%
output_matrix_gain[which(output_matrix_gain > 100)] <- 100


common_samples <- intersect(row.names(output_matrix_gain), row.names(output_matrix_loss))

output_matrix_loss <- output_matrix_loss[common_samples,]
output_matrix_gain <- output_matrix_gain[common_samples,]

#CNA calcs (cutoff > 0% of chromosome arm)

loss_mtx <- output_matrix_loss
gain_mtx <- output_matrix_gain
all_mtx <- loss_mtx + gain_mtx

cell_losses <- sapply(1:nrow(loss_mtx),function(row){sum(loss_mtx[row,] > 0)})
cell_losses_plot_data <- data.frame(
  count = cell_losses
  ,sample = "CCLE"
  ,type = "losses"
)

cell_gains <- sapply(1:nrow(gain_mtx),function(row){sum(gain_mtx[row,] > 0)})
cell_gains_plot_data <- data.frame(
  count = cell_gains
  ,sample = "CCLE"
  ,type = "gains"
)

cell_all <- sapply(1:nrow(all_mtx),function(row){sum(all_mtx[row,] > 0)})
cell_all_plot_data <- data.frame(
  count = cell_all
  ,sample = "CCLE"
  ,type = "all"
)


cell_plot_data <- rbind(cell_losses_plot_data, cell_gains_plot_data, cell_all_plot_data)

#aneuploidy calcs (cutoff > value of chromosome arm)

aneuploid_losses <- sapply(1:nrow(loss_mtx),function(row){sum(loss_mtx[row,] > percent_of_arm_threshold)})
aneuploid_losses_plot_data <- data.frame(
  count = aneuploid_losses
  ,sample = "CCLE"
  ,type = "losses"
)

aneuploid_gains <- sapply(1:nrow(gain_mtx),function(row){sum(gain_mtx[row,] > percent_of_arm_threshold)})
aneuploid_gains_plot_data <- data.frame(
  count = aneuploid_gains
  ,sample = "CCLE"
  ,type = "gains"
)

aneuploid_all <- sapply(1:nrow(all_mtx),function(row){sum(all_mtx[row,] > percent_of_arm_threshold)})
aneuploid_all_plot_data <- data.frame(
  count = aneuploid_all
  ,sample = "CCLE"
  ,type = "all"
)


aneuploid_plot_data <- rbind(aneuploid_losses_plot_data, aneuploid_gains_plot_data, aneuploid_all_plot_data)

output_matrix_gains_thresholded <- output_matrix_gain
output_matrix_gains_thresholded[which(output_matrix_gains_thresholded < percent_of_arm_threshold)] <- 0 
output_matrix_gains_thresholded[which(output_matrix_gains_thresholded >= percent_of_arm_threshold)] <- 1

output_matrix_losses_thresholded <- output_matrix_loss
output_matrix_losses_thresholded[which(output_matrix_losses_thresholded < percent_of_arm_threshold)] <- 0 
output_matrix_losses_thresholded[which(output_matrix_losses_thresholded >= percent_of_arm_threshold)] <- 1

output_matrix_signed_thresholded <- output_matrix_gains_thresholded - output_matrix_losses_thresholded

#remove mixed chromosomes, set as not aneuploid. Only needed for some low thresholds chosen by user
output_matrix_gains_thresholded[
which( (output_matrix_gains_thresholded == 1) & (output_matrix_losses_thresholded==1))
] <- 0

output_matrix_losses_thresholded[
  which( (output_matrix_gains_thresholded == 1) & (output_matrix_losses_thresholded==1))
] <- 0


output_matrix_signed_thresholded <- output_matrix_signed_thresholded[,order(colnames(output_matrix_signed_thresholded))]

write.table(output_matrix_signed_thresholded,output_arm_CNAs,col.names = NA, row.names=TRUE, sep="\t" )



