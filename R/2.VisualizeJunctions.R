# Author:       Job van Riet
# Date:         09-01-23
# Function:     Visualize splice-junctions detected within mCRPC samples.

library(plyr)
library(dplyr)
library(patchwork)

source('R/misc_Themes.R')

# Load data. ----

DR71.ARvs <- readRDS('/mnt/share1/repository/HMF/DR71/Dec2021/RData/DR71.ARvs.Rds')
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/DR71.RNASeq.RData')

# Retrieve AR gene-expression. ----

counts.AR <- DESeq2::plotCounts(DR71.RNASeq$DESeq2.withoutSequencingBatch, gene = 'ENSG00000169083', intgroup = 'sequencingBatch', returnData = T, normalized = T)
counts.AR$sample <- DR71.RNASeq$DESeq2.withoutSequencingBatch$sample


# Determine expressed ARVs. ----

AR.PSI <- reshape2::melt(DR71.ARvs$AR.PSI) %>%
    dplyr::mutate(value = ifelse(is.nan(value), 0, value)) %>%
    dplyr::filter(grepl('PSI', variable)) %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(
        median = median(value, na.rm = T),
        maxScore = max(value),
        nSamplesWithExpression = sum(value >= 0.01),
        variantId = gsub('3.7kbUTR', '<sup>3.7kbUTR</sup>', gsub('v', '-V', gsub('PSI.', '', variable))),
        iqr = sprintf('%s (%s-%s)', round(median * 100, 2), round(summary(value)[1] * 100, 2), round(summary(value)[4] * 100, 2)),
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(maxScore > 0) %>% 
    dplyr::filter(!grepl('UTR', variantId))


# Visualize relative expression of variants. -----

plot.N <- AR.PSI %>% dplyr::distinct(variantId, median, nSamplesWithExpression, variable) %>% 
    ggplot2::ggplot(., ggplot2::aes(x = reorder(variantId, -median), y = nSamplesWithExpression, fill = variable, label  = nSamplesWithExpression)) +
    ggplot2::geom_bar(stat = 'identity', color = 'black', lwd = .5, width = .8) +
    ggplot2::geom_text(nudge_y = +10, size = 3) +
    ggplot2::scale_y_continuous(breaks = c(seq(0, 200, 25)), expand = c(0,.003), limits = c(0, 200)) +
    ggplot2::labs(x = NULL, y = 'No. samples<br>(Rel. expression ≥ 1%)') +
    ggplot2::scale_fill_grey(guide = 'none') +
    theme_Job +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank())

plot.PSI <- ggplot2::ggplot(AR.PSI, ggplot2::aes(x = reorder(variantId, -median), y = value, fill = variable)) +
    ggplot2::geom_hline(yintercept = 0.01, color = 'red', lty = 'dotted') +
    gghalves::geom_half_boxplot(outlier.shape = NA) +
    gghalves::geom_half_point_panel(shape = 21, position = ggbeeswarm::position_quasirandom(method = 'tukey', width = .15)) +
    ggplot2::scale_y_sqrt(labels = scales::percent_format(), breaks = c(0, .01, seq(0, .2, .025)), expand = c(0,.003), limits = c(0, .2)) +
    ggplot2::labs(x = 'AR splice-variants<br>(CPCT-02; <i>n</i> = 278)', y = 'Rel. expression of variant-specific SJ<br>(SJ<sub><i>x</i></sub> / SJ<sub>E1E2</sub>)') +
    ggplot2::scale_fill_grey(guide = 'none') +
    theme_Job +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank())

plot.N + plot.PSI + patchwork::plot_layout(ncol = 1, heights = c(.5, 1))


# Heatmap of AR-Vs per sample. ----

dataPSI <- AR.PSI %>% 
    dplyr::mutate(value = ifelse(value >= 0.01, T, F)) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::mutate(totalN = sum(value)) %>% 
    dplyr::ungroup()


heatData <- reshape2::dcast(data = dataPSI, formula = sample ~ variantId, value.var = 'value')
rownames(heatData) <- heatData$sample; heatData$sample <- NULL

## Perform memosort and order genes/samples. ----

heatData <- dataPSI %>% 
    dplyr::inner_join(counts.AR %>% dplyr::distinct(sample, count)) %>% 
    dplyr::mutate(
        hasARv = ifelse(totalN > 0, 'Yes', 'No'),
        totalN = ifelse(totalN >=5, '≥5', totalN),
        totalN = factor(totalN, levels = c(0:4, '≥5')),
        variantId = factor(variantId, levels = rev(rownames(R2CPCT::memoSort(as.data.frame(t(heatData)))))),
        sample = factor(sample, colnames(R2CPCT::memoSort(as.data.frame(t(heatData)))))
    )

plot.heat <- ggplot2::ggplot(heatData, ggplot2::aes(x=sample, y = variantId, fill = value)) +
    ggplot2::geom_tile(color = "grey95", lwd = 0.2, linetype = 1) +
    ggplot2::labs(x = 'CPCT-02 (<i>n</i> = 278)', y = 'Observed AR-Vs<br>(Rel. expression ≥ 1%)') +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_manual(values = c('white', 'black'), guide = 'none') +
    theme_Job2

plot.ar <- ggplot2::ggplot(heatData %>% dplyr::distinct(sample, count), ggplot2::aes(x=sample, xend = sample, yend = 0, y = count)) +
    ggplot2::geom_segment(color = 'grey50') +
    ggplot2::geom_point(size = .5) +
    ggplot2::labs(x = NULL, y = 'Expression - <i>AR</i><br>(Normalized)') +
    ggplot2::scale_y_sqrt(breaks = c(0, 10000, 50000, 100000, 250000, 500000, 1000000), limits = c(0, 1000000), expand = c(0,0)) +
    theme_Job2 + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = 'grey75', linetype = '12'))

plot.ar + plot.heat + patchwork::plot_layout(ncol = 1, heights = c(1, .66))


# Visualize ARwt vs. no. detected ARVs ----

stat.test <- heatData %>% 
    dplyr::distinct(count, totalN, hasARv) %>% 
    rstatix::pairwise_wilcox_test(count ~ totalN)

stat.test2 <- heatData %>% 
    dplyr::distinct(count, totalN, hasARv) %>% 
    rstatix::pairwise_wilcox_test(count ~ hasARv)

plot.totalN <- heatData %>% 
    dplyr::distinct(count, totalN, sample) %>% 
    ggplot2::ggplot(., ggplot2::aes(x = totalN, y = count, group = totalN, fill = totalN)) + 
    gghalves::geom_half_boxplot(outlier.shape = NA) +
    gghalves::geom_half_point_panel(shape = 21, position = ggbeeswarm::position_quasirandom(method = 'tukey', width = .1)) +
    ggplot2::scale_y_sqrt(breaks = c(0, 10000, 50000, 100000, 250000, 500000, 1000000), limits = c(0, 1000000), expand = c(0,002)) +
    ggplot2::scale_fill_grey(guide = 'none') +
    ggplot2::labs(x = 'No. of ARVs<br>(Rel. expression ≥ 1%)', y = 'Expression - <i>AR</i> (Normalized)') +
    ggpubr::geom_bracket(data = stat.test %>% dplyr::filter(p.adj.signif != 'ns'), inherit.aes = F, mapping = ggplot2::aes(xmin = group1, xmax = group2), y.position = 800, step.increase = 0.04, tip.length = 0.01, family = 'Nimbus Sans') +
    theme_Job

plot.anyV <- heatData %>% 
    dplyr::distinct(count, hasARv, sample) %>% 
    ggplot2::ggplot(., ggplot2::aes(x = hasARv, y = count, group = hasARv, fill = hasARv)) + 
    gghalves::geom_half_boxplot(outlier.shape = NA) +
    gghalves::geom_half_point_panel(shape = 21, position = ggbeeswarm::position_quasirandom(method = 'tukey', width = .1)) +
    ggplot2::scale_y_sqrt(breaks = c(0, 10000, 50000, 100000, 250000, 500000, 1000000), limits = c(0, 1000000), expand = c(0,002)) +
    ggplot2::scale_fill_grey(guide = 'none') +
    ggplot2::labs(x = 'Presence of any ARv<br>(Rel. expression ≥ 1%)', y = 'Expression - <i>AR</i> (Normalized)') +
    ggpubr::geom_bracket(data = stat.test2 %>% dplyr::filter(p.adj.signif != 'ns'), inherit.aes = F, mapping = ggplot2::aes(xmin = group1, xmax = group2), y.position = 800, step.increase = 0.02, tip.length = 0.01, family = 'Nimbus Sans') +
    theme_Job

plot.anyV + plot.totalN + patchwork::plot_layout(ncol = 2, widths = c(1, 2))
