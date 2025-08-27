# rm(list = ls())
library(magrittr)

# tf----
load("B-network.rds")

edges <- mrna_tf %>% 
  data.frame(check.names = F) %>% 
  `names<-`(c("from", "to")) %>% 
  dplyr::distinct()
nodes <- tf_attr %>% 
  data.frame(check.names = F) %>% 
  `names<-`(c("id", "type"))
nodes$type <- factor(nodes$type, levels = c("mRNA", "TF"))

library(tidygraph)
net.tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

library(ggraph)
set.seed(123)
ggraph(net.tidy, layout = "graphopt") + 
  geom_node_point(aes(color = type, fill = type, shape = type, size = type)) + # 点信息
  geom_edge_link(alpha = 1, width = 0.1) +  # 边信息
  # scale_edge_width(range = c(0.2, 2)) + # 控制粗细
  geom_node_text(aes(label = id), repel = TRUE, check_overlap = T, size = 5) + # 增加节点的标签，reple避免节点重叠
  # labs(edge_width  = "score") + # 图例标签
  scale_size_manual(values = c(6, 3)) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("#008B45", "#3B4992")) + 
  scale_color_manual(values = c("#008B45", "#3B4992")) + 
  theme_graph() +
  coord_cartesian(clip = "off")

ggsave("tf.pdf", width = 8, height = 8, device = cairo_pdf)

# mirna----
edges <- mRNA_miRNA %>% 
  data.frame(check.names = F) %>% 
  `names<-`(c("from", "to")) %>% 
  dplyr::distinct()
nodes <- miRNA_attr %>% 
  data.frame(check.names = F) %>% 
  `names<-`(c("id", "type"))
nodes$type <- factor(nodes$type, levels = c("mRNA", "miRNA"))

library(tidygraph)
net.tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

library(ggraph)
set.seed(123)
ggraph(net.tidy, layout = "graphopt") + 
  geom_node_point(aes(color = type, fill = type, shape = type, size = type)) + # 点信息
  geom_edge_link(alpha = 1, width = 0.1) +  # 边信息
  # scale_edge_width(range = c(0.2, 2)) + # 控制粗细
  geom_node_text(aes(label = id), size = 5) + # 增加节点的标签，reple避免节点重叠
  # labs(edge_width  = "score") + # 图例标签
  scale_size_manual(values = c(6, 3)) +
  scale_shape_manual(values = c(21, 23)) +
  scale_fill_manual(values = c("#008B45", "#EE0000")) + 
  scale_color_manual(values = c("#008B45", "#EE0000")) + 
  theme_graph() +
  coord_cartesian(clip = "off")
ggsave("mirna.pdf", width = 8, height = 8, device = cairo_pdf)
