# Author:       Job van Riet
# Date:         04-01-23
# Function:     Detect and quantify splice-junctions within the gene-body of AR within mCRPC samples.

library(plyr)
library(dplyr)

source('R/misc_Functions.R')

# Import meta-data of DR-071 samples ----

samples.meta <- readxl::read_xlsx('Misc/SupplTable1_OverviewData.xlsx', sheet = 'Samples (HMF)')

# Retrieve BAM files ----

# Select correct samples.
files.BAM <- data.frame(file_bam = list.files('/mnt/share1/repository/HMF/DR71/Dec2021/dataHMF/RNASeq/BAM/', pattern = 'Aligned.out.sorted.markDup.bam$', full.names = T)) %>%
    dplyr::mutate(
        sample = gsub('_.*', '', basename(file_bam)),
        paired_end = T
    ) %>%
    dplyr::inner_join(samples.meta, by = 'sample')


# Detect splice-junctions surrounding AR ----

retrieveJunctions <- function(f, param, strandModus){
    
    x <- Rsamtools::BamFile(f, asMates = T)
    
    # Count junction reads.
    ga <- GenomicAlignments::readGAlignmentPairs(x, param = param, strandMode = strandModus)
    GenomeInfoDb::seqlevelsStyle(ga) <- 'UCSC'
    ga <- GenomeInfoDb::keepSeqlevels(ga, value = GenomeInfoDb::seqlevelsInUse(ga), pruning.mode = 'coarse')
    
    # Summarize reads per junction and retrieve the intron motif.
    ga <- GenomicAlignments::summarizeJunctions(ga, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
    
    # Add sample information.
    ga$sample <- gsub('_.*', '', basename(f))
    
    # Return data.
    return(ga)
}

# Only count primary, QC-passed, properly-paired and non-duplicated reads.
param <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges('X', IRanges::IRanges(66763863,66951462)), flag = Rsamtools::scanBamFlag(isPaired = T, isProperPair = T, isSecondaryAlignment = F, isNotPassingQualityControls = F, isDuplicate = F))

AR.junctions <- unlist(GenomicRanges::GRangesList(pbapply::pbapply(files.BAM, 1, function(x){
    
    # Determine correct strand mode based on sequencing platform/kit.
    strandModus <- ifelse(x['SequencingPlatform'] == "Illumina NovaSeq (~250bp)", 1, 2)
    return(retrieveJunctions(x['file_bam'], param, strandModus))
    
}, cl = 30)))

# Annotate the junctions.
# These should be 1nt away from the known (cryptic) exons. Otherwise, there's a junction to a novel exon and will map to the 'Unknown' element.
anno <- IRanges::follow(AR.junctions, geneStructure.AR)
anno <- ifelse(is.na(anno), length(geneStructure.AR), anno)
anno <- ifelse(IRanges::start(AR.junctions) - IRanges::end(geneStructure.AR[anno]) == 1, anno, length(geneStructure.AR))
AR.junctions$Exon.Donor <- geneStructure.AR[anno]$Name

anno <- IRanges::precede(AR.junctions, geneStructure.AR)
anno <- ifelse(is.na(anno), length(geneStructure.AR), anno)
anno <- ifelse(IRanges::end(AR.junctions) - IRanges::start(geneStructure.AR[anno]) == -1, anno, length(geneStructure.AR))
AR.junctions$Exon.Acceptor <- geneStructure.AR[anno]$Name


# Filter on relevant splice-junctions ----

# Convert to tibble.
AR.junctions.df <- tibble::as_tibble(base::data.frame(AR.junctions)) %>%
    
    # Min. spliced-reads in sample.
    dplyr::filter(score >= 5) %>% 
    
    # Per junction, count the number of samples (with enough reads).
    dplyr::group_by(Exon.Donor, Exon.Acceptor, start, end, width, intron_motif) %>%
    dplyr::mutate(totalSamples = dplyr::n_distinct(sample), meanJunctionReadsAllSamples = mean(score)) %>%
    dplyr::ungroup() %>%
    
    # Only select junctions found in multiple samples.
    dplyr::filter(totalSamples >= 10) %>% 
    dplyr::mutate(sample = factor(sample, levels = unique(AR.junctions$sample))) %>% 
    tidyr::complete(sample, fill = list(0))

# Determine unknown junctions.
# AR.junctions.df %>% dplyr::filter(Exon.Donor == 'Unknown' | Exon.Acceptor == 'Unknown') %>% View


# Determine percent spliced in (PSI) compared to ARwt ----

# Determine ratio of AR-Vs per sample compared to their counterpart wild-type upstream junction.
AR.PSI <-  base::do.call(base::rbind, lapply(base::split(AR.junctions.df, AR.junctions.df$sample), function(x){
    
    isEmpty <- function(x) {
        return(length(x) == 0)
    }
    
    getReads <- function(x, ExonA, ExonB){
        z <- (x %>% dplyr::filter(Exon.Donor == ExonA, Exon.Acceptor == ExonB))$score
        z <- ifelse(isEmpty(z), 0, z)
        return(z)
    }
    
    z <- tibble::tibble(
        
        # Junction-counts for wt-junctions of AR.
        counts.AR.wt.1_2 = getReads(x, 'Exon 1', 'Exon 2'),
        counts.AR.wt.2_3 = getReads(x, 'Exon 2', 'Exon 3'),
        counts.AR.wt.3_4 = getReads(x, 'Exon 3', 'Exon 4'),
        counts.AR.wt.4_5 = getReads(x, 'Exon 4', 'Exon 5'),
        counts.AR.wt.6_7 = getReads(x, 'Exon 6', 'Exon 7'),
        counts.AR.wt.7_8 = getReads(x, 'Exon 7', 'Exon 8'),
        
        # Junction-counts for the AR-Vs.
        counts.ARv1 = max(getReads(x, 'Exon 3', 'CE1'), getReads(x, 'Exon 3', 'CE1\'')),
        counts.ARv2 = getReads(x, 'Exon 3', 'Exon 3'),
        counts.ARv3 = getReads(x, 'Exon 2', 'CE4'),
        counts.ARv4 = getReads(x, 'Exon 3', 'CE4'),
        counts.ARv5 = getReads(x, 'Exon 3', 'CE2'),
        counts.ARv6 = getReads(x, 'Exon 3', 'CE2*'),
        counts.ARv7 = getReads(x, 'Exon 3', 'CE3'),
        counts.ARv8 = getReads(x, 'Exon 3', 'CTE1'),
        counts.ARv9 = getReads(x, 'Exon 3', 'CE5'),
        counts.ARv10 = getReads(x, 'Exon 3', 'CTE2'),
        counts.ARv12 = getReads(x, 'Exon 4', 'Exon 8*'),
        counts.ARv13 = getReads(x, 'Exon 6', 'Exon 8*'),
        counts.ARv14 = getReads(x, 'Exon 7', 'Exon 8*'),
        counts.AR8 = getReads(x, 'Exon 1', 'CTE1'),
        counts.AR23 = getReads(x, 'Exon 2', 'Exon 3 - Extension69'),
        counts.AR45 = getReads(x, 'Exon 1b', 'Exon 2'),
        
        # UTR variants.
        counts.3.7kbUTR = getReads(x, 'UTR Shortening - 3.7kb-UTR', '3\' UTR*'),
        counts.0.9kbUTR = getReads(x, 'UTR Shortening - 0.9kb-UTR', '3\' UTR*'),

        # Derive the PSI by comparing against E1-E2.
        PSI.ARv1 = counts.ARv1 / counts.AR.wt.1_2,
        PSI.ARv2 = counts.ARv2 / counts.AR.wt.1_2,
        PSI.ARv3 = counts.ARv3 / counts.AR.wt.1_2,
        PSI.ARv4 = counts.ARv4 / counts.AR.wt.1_2,
        PSI.ARv5 = counts.ARv5 / counts.AR.wt.1_2,
        PSI.ARv6 = counts.ARv5 / counts.AR.wt.1_2,
        PSI.ARv7 = counts.ARv7 / counts.AR.wt.1_2,
        PSI.ARv8 = counts.ARv8 / counts.AR.wt.1_2,
        PSI.ARv9 = counts.ARv9 / counts.AR.wt.1_2,
        PSI.ARv10 = counts.ARv10 / counts.AR.wt.1_2,
        PSI.ARv12 = counts.ARv12 / counts.AR.wt.1_2,
        PSI.ARv13 = counts.ARv13 / counts.AR.wt.1_2,
        PSI.ARv14 = counts.ARv14 / counts.AR.wt.1_2,
        PSI.AR8 = counts.AR8 / counts.AR.wt.1_2,
        PSI.AR23 = counts.AR23 / counts.AR.wt.1_2,
        PSI.AR45 = counts.AR45 / counts.AR.wt.1_2,
        PSI.AR3.7kbUTR = counts.3.7kbUTR / counts.AR.wt.1_2,
        PSI.AR0.9kbUTR = counts.0.9kbUTR / counts.AR.wt.1_2,
        
        sample = unique(x$sample)
    )
    
    return(z)
}))

# Save data ----

DR71.ARvs <- list(AR.junctions = AR.junctions, AR.junctions.df = AR.junctions.df, AR.PSI = AR.PSI)
saveRDS(DR71.ARvs, file = '/mnt/share1/repository/HMF/DR71/Dec2021/RData/DR71.ARvs.Rds')
