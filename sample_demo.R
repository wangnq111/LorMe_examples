library(LorMe)
library(magrittr)
library(ggplot2)

data("Two_group")

LorMe_defaults()
LorMe_options(
  beta=list(diagram ="ellipse"),
  global = list(Analysis_level = "Genus"),
  sub_net = list(threshold = 0.6,n=3,method="spearman",reads=T),
  all_net = list(threshold = 0.75,n=2,method="pearson",reads=T)
)
log <- system.time({Results=LorMe_pipeline(Two_group)}) #,step = "all_net"


((Results$alpha_results$plotlist$Plotobj_Shannon$Barplot|
Results$beta_results$PCoA_Plot|
Results$composition_results$alluvialplot)/
((Results$Deseq_volcano$Mean_FC|
Results$Deseq_manhattan$manhattan)+plot_layout(widths = c(0.75,2))))+plot_layout(guides = "collect")
ggsave("./example.pdf",width = 8,height = 4)
#Results$sub_network_results$Treatment_sub_network %>% network_visual()

pdf("./example_network.pdf",width = 8,height = 3)
Results$combine_network_results %>% network_visual(major_num = 10)
dev.off()

pdf("./example_network_2.pdf",width = 8,height = 3)
Results$combine_network_results %>% network_visual(mode = "major_tax",
                                                   taxlevel = "Phylum",select_tax =names(Results$composition_results$filled_color)[1:10],palette = Results$composition_results$filled_color)
dev.off()


