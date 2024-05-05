d2 = read.csv("../../data/sample_counts.csv")
comet_cols = pals::brewer.reds(12)[c(2,5,8,11)]
other_cols = pals::kelly(6)[c(3,6)]
pdf("sample_count.pdf", width=10, height=6)
ggplot(d2, aes(paste(disease,"\n",steroid), count, fill=paper )) + geom_bar(stat = "identity", position = "stack") + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust = 1)) + facet_grid(. ~ paste(tissue,"\n",modality)) + scale_fill_manual(values = c(comet_cols[c(1,4,2)], other_cols[1], comet_cols[3], other_cols[2])) + ylab("Sample count")
dev.off()

