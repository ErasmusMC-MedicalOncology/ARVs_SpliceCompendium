
# Enhancer & gene-body structure of AR (GRCh37) ----

# Additional (cryptic) exons were denoted based on the presence of >=10 splice-junctions in >=10 mCRPC samples in the CPCT-02 cohort.

geneStructure.AR <- unlist(GenomicRanges::GRangesList(
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66132562, 66134560), Type = 'Enhancer', Name = 'Promoter/Enhancer #1', Comment = 'van Dessel / van Riet et al.'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66763221, 66771842), Type = 'Enhancer', Name = 'Promoter/Enhancer #2', Comment = 'GH0XJ067543'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66763863, 66764149), Type = 'UTR', Name = '5\' UTRs'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66763863, 66764464), Type = 'UTR', Name = '5\' UTR'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66764465, 66767866), Type = 'Exon', Name = 'Exon 1l'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66764465, 66766604), Type = 'Exon', Name = 'Exon 1'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66788633, 66788864), Type = 'Exon', Name = 'Exon 1b', Comment = 'Alternative Exon 1'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66788678, 66788678+150), Type = 'Exon', Name = 'CE6a', Comment = 'Cryptic Exon 6 - Job'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66788681, 66788681+150), Type = 'Exon', Name = 'CE6b', Comment = 'Cryptic Exon 6 - Job'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66861773-150, 66861773), Type = 'Exon', Name = 'CEu2', Comment = 'Unknown Cryptic Exon upstream of E2'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66863098, 66863249), Type = 'Exon', Name = 'Exon 2'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66900562, 66900845), Type = 'Exon', Name = 'CE4', Comment = 'Cryptic Exon'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66905783, 66905851), Type = 'Exon', Name = 'Exon 3 - Extension69', Comment = 'Extension of Exon 3 (69nt)'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66905852, 66905968), Type = 'Exon', Name = 'Exon 3'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66906802, 66906802 + 150), Type = 'Exon', Name = 'CEu1', Comment = 'Cryptic Exon (Unknown)'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66909398, 66909398 + 150), Type = 'Exon', Name = 'CE1',  Comment = 'Cryptic Exon'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66909401, 66909401 + 150), Type = 'Exon', Name = 'CE1\'', Comment = 'Cryptic Exon'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66910123, 66910123 + 150), Type = 'Exon', Name = 'CTE1', Comment = 'Cryptic Stop Sequence'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66912029, 66912029 + 150), Type = 'Exon', Name = 'CE2', Comment = 'Cryptic Exon'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66912110, 66912110 + 150), Type = 'Exon', Name = 'CE2*', Comment = 'Cryptic Exon'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66913412, 66913412 + 150), Type = 'Exon', Name = 'CE5', Comment = 'Cryptic Exon'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66914515, 66915918), Type = 'Exon', Name = 'CE3', Comment = 'Cryptic Exon'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66914721, 66914721 + 150), Type = 'Exon', Name = 'CTE2', Comment = 'Cryptic Stop Sequence'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66931244, 66931531), Type = 'Exon', Name = 'Exon 4'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66937320, 66937464), Type = 'Exon', Name = 'Exon 5'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66941675, 66941805), Type = 'Exon', Name = 'Exon 6'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66942669, 66942826), Type = 'Exon', Name = 'Exon 7'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66943528, 66943683), Type = 'Exon', Name = 'Exon 8'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66943608, 66943683), Type = 'Exon', Name = 'Exon 8*'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66945478-150, 66945478), Type = 'Exon', Name = 'UTR Shortening - 3.7kb-UTR'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66947623 - 150, 66947623), Type = 'Exon', Name = 'UTR Shortening - 5.9kb-UTR'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66948645, 66948645 + 150), Type = 'Exon', Name = 'UTR Shortening - 3.6kb-UTR'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66948516, 66950461), Type = 'Exon', Name = '3\' UTR*'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(66950462, 66951462), Type = 'UTR', Name = '3\' UTR'),
    GenomicRanges::GRanges('chrX', IRanges::IRanges(1, 2), Type = 'Other', Name = 'Unknown')
))

GenomicRanges::strand(geneStructure.AR) <- '+'
geneStructure.AR$gene_name <- 'AR'