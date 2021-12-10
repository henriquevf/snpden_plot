library(ggplot2)
library(gdata)

snpden <- read.table('neofelis_nebulosa.g.hetSites.snpden', header = T)

target <- c('HiC_scaffold_1','HiC_scaffold_4','HiC_scaffold_9',
            'HiC_scaffold_3','HiC_scaffold_6','HiC_scaffold_7',
            'HiC_scaffold_8','HiC_scaffold_2','HiC_scaffold_5',
            'HiC_scaffold_11','HiC_scaffold_14','HiC_scaffold_13','HiC_scaffold_12',
            'HiC_scaffold_18','HiC_scaffold_17','HiC_scaffold_19','HiC_scaffold_16',
            'HiC_scaffold_15','HiC_scaffold_10')

chr <-c('A1','A2','A3','B1','B2','B3','B4','C1','C2','D1','D2','D3','D4',
        'E1','E2','E3','F1','F2','X')

chr.len <-c('0Mbp','50Mbp','100Mbp','150Mbp','200Mbp','250Mbp','300Mbp')

snpden$CHROM <- reorder.factor(snpden$CHROM, new.order = target)
snpden <-subset(snpden, snpden$CHROM!='NA')

snpden$groups = cut(snpden$VARIANTS.KB, c(0,0.05,0.1,0.15,0.20,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,
                                          3,3.25,3.5,3.75,4,4.25,4.5,4.75,5))

levels(snpden$CHROM) <-c('A1','A2','A3','B1','B2','B3','B4','C1','C2','D1','D2','D3','D4',
                         'E1','E2','E3','F1','F2','X')

title<-expression(paste(italic("Neofelis nebulosa"), " heterozygous SNP densities"))

snpden_plot <- ggplot(data = snpden, aes(x=BIN_START, y=1)) + 
  geom_tile(aes(fill=groups)) +
  facet_grid(CHROM ~ ., switch='y') +
  xlab('Chromosome Length') + ylab('Scaffold Number') + theme_minimal() +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 180)) +
  scale_fill_manual(values = rainbow(21, start = .4, end = .1, rev = FALSE,),
                    name = "Variants/kb", 
                    labels = c("<0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25",
                               "0.25-0.50","0.50-0.75","0.75-1.0","1.0-1.25","1.25-1.5",
                               "1.5-1.75","1.75-2.0","2.0-2.25","2.25-2.5","2.5-2.75","2.75-3.0",
                               "3.0-3.25","3.25-3.5","3.5-3.75","3.75-4.0","4.0-4.25","4.25-4.5",
                               "4.5-4.75","4.75-5")) +
  scale_x_continuous(name='Chromosome length', labels = c('0Mbp','50Mbp','100Mbp',
                                                          '150Mbp','200Mbp','250Mbp'))


ggsave('Nnebulosa.1Mb.snpden.v2.svg', plot = snpden_plot, device = 'svg',
       dpi = 600, units = c('cm'), width = 20, height = 18)

ggsave('Nnebulosa.1Mb.snpden.v2.pdf', plot = snpden_plot, device = 'pdf',
       dpi = 600, units = c('cm'), width = 20, height = 18)
