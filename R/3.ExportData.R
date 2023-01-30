# Author:       Job van Riet
# Date:         19-01-23
# Function:     Export data into Suppl. Table.

library(plyr)
library(dplyr)

DR71.ARvs <- readRDS('/mnt/share1/repository/HMF/DR71/Dec2021/RData/DR71.ARvs.Rds')
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/DR71.RNASeq.RData')

# Retrieve AR gene-expression. ----

counts.AR <- DESeq2::plotCounts(DR71.RNASeq$DESeq2.withoutSequencingBatch, gene = 'ENSG00000169083', intgroup = 'sequencingBatch', returnData = T, normalized = T)
counts.AR$sample <- DR71.RNASeq$DESeq2.withoutSequencingBatch$sample

# Open workbook.
wb <- openxlsx::createWorkbook(creator = 'Job van Riet', subject = 'AR-Vs', title = 'AR-Vs', category = 'CPCT-02')

# Expression of AR-Vs - CPCT-02 ----
openxlsx::addWorksheet(wb, sheetName = 'mCRPC (CPCT-02)')

data.ARv <- reshape2::melt(DR71.ARvs$AR.PSI) %>%
    dplyr::mutate(value = ifelse(is.nan(value), 0, value)) %>% 
    reshape2::dcast(data = ., formula = sample ~ variable, value.var = 'value') %>% 
    dplyr::inner_join(data.frame(SummarizedExperiment::colData(DR71.RNASeq$DESeq2.withoutSequencingBatch)) %>% dplyr::distinct(sample, hmfSampleId)) %>% 
    dplyr::inner_join(counts.AR) %>% 
    dplyr::mutate(sample = hmfSampleId, hmfSampleId = NULL, sequencingBatch = NULL) %>% 
    dplyr::select(sample, countAR = count, dplyr::everything())


openxlsx::writeDataTable(wb, sheet = 'mCRPC (CPCT-02)', x = data.ARv)
openxlsx::saveWorkbook(wb, file = 'SupplTable1_OverviewData.xlsx', overwrite = T)
