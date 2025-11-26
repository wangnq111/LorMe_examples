log <- system.time({
library(LorMe)
library(magrittr)
library(patchwork)
library(ggplot2)


otu=read.table("./otu_table.txt",header = T,sep="\t")
design=read.table("./design.txt",header = T,sep="\t")



my_color=color_scheme("Plan3") %>% .[c(1,6)]
names(my_color)=c("Sensitive","Resistance")
Lorme_object=tax_summary(groupfile = design,inputtable =otu[,c(2:33)],reads = T,taxonomytable =otu[,c(1,34)],outputtax = "standard")
Lorme_object_config=object_config(taxobj =Lorme_object,treat_location = 2,treat_col =  my_color,treat_order = c("Sensitive","Resistance"))



LorMe_options(
  global = list(Analysis_level = "Genus"),
  comp=list(palette="Paired"),
  deseq=list(control_name="Sensitive"),
  sub_net = list(threshold = 0.75,n=8),
  all_net = list(threshold = 0.75,n=8)
)

results_all=LorMe_pipeline(Lorme_object_config, step = "all")

phylum_col=results_all$composition_results$filled_color
  


####
a=results_all$alpha_results$plotlist$Plotobj_Shannon$Violinplot
b=results_all$beta_results$PCoA_Plot
c=results_all$composition_results$barplot
line1=(a|b|c)+plot_layout(widths = c(1,1,2),guides = "collect")
line1
ggsave("./line1.pdf",width = 13,height = 3.3)

d=results_all$indic_volcano$FC_FDR
e=results_all$indic_manhattan$manhattan_circle
(d|e)+plot_layout(,guides = "collect")
ggsave("./line2_left.pdf",width = 6.5,height = 3.3)

pdf("./line2.pdf",width = 13,height = 3.3)
par(mfrow = c(1, 2), mar = c(0, 1, 1, 0), pty = "m")
s_net=network_visual(network_obj = results_all$sub_network_results$Sensitive_sub_network,mode = "major_tax",major_num = 5,taxlevel = "Phylum",select_tax =names(phylum_col)[1:10],palette = phylum_col )
r_net=network_visual(network_obj = results_all$sub_network_results$Resistance_sub_network,mode = "major_tax",major_num = 5,taxlevel = "Phylum",select_tax =names(phylum_col)[1:10],palette = phylum_col )
dev.off()

network_all=results_all$combine_network_results
pdf("./line3_left.pdf",width = 13,height = 3.3)
meta_net=network_withdiff(network_obj =network_all,
                          diff_frame =results_all$Deseq_results,
                          aes_col =   Lorme_object_config$configuration$treat_col)
dev.off()

Module_abundance=Module_abundance(network_obj = network_all,No.module = c(2,3))
Module_composition=Module_composition(network_obj =network_all,No.module = c(2,3),taxlevel = "Phylum",mode = "select",palette = phylum_col,select_tax = names(phylum_col)[1:10] )

h1=Module_abundance$plotlist$Plotobj_Module2$Violinplot+guides(fill="none")
h2=Module_composition$Module2$Pie_plot_Module2+guides(fill="none")
h3=Module_abundance$plotlist$Plotobj_Module3$Violinplot+guides(fill="none")
h4=Module_composition$Module3$Pie_plot_Module3+guides(fill="none")
((h1|h2)/(h3|h4))
ggsave("./line3_right.pdf",width = 3,height = 3.3)

})
