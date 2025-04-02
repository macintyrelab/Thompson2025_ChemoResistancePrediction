# Clean environment
freshr::freshr()

# Libraries
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(this.path)
library(dplyr)

# Paths
pancancer_dir <- dirname(this.path())
base_dir <- dirname(dirname(pancancer_dir))
figs_dir <- file.path(base_dir, 'Figures')
supp_figs_dir <- file.path(figs_dir, 'Supplementary')
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_figs_dir, showWarnings = FALSE, recursive = TRUE)

# Render graph
ndf <- data.frame(id=1:17, type=c('All','HMF','TCGA','TCGA','Breast','Ovarian',
                                  'All','All','All','All','All','All','All',
                                  'All','All', 'All','TCGA'),
                  label=c("Filter to cancer type",
                          "Keep only treatment lines\nimmediately post-biopsy",
                          "Remove treatment lines without TTF values",
                          "Remove Stage I patients",
                          "Remove patients treated\nwith Anti-HER2s",
                          "Remove 1st line treatments",
                          "Predict chemotherapy resistance",
                          "Identify experimental arm treatment lines",
                          "Identify control arm treatment lines",
                          "Remove treatment lines with insufficient cycles",
                          "Remove treatment lines with insufficient cycles",
                          "Remove overlapping patients from control arm",
                          "Identify most common SoC", "Keep only most common SoC",
                          "Keep only 1 treatment line per patient",
                          "Keep only 1 treatment line per patient",
                          "Exclude patients with long followup periods"))
edf <- data.frame(id=1:21,
                  from=c(1, 1, 2, 2, 3, 4, 4, 5, 4, 6, 7, 7, 8, 9, 11, 12, 13, 10, 14, 15, 16),
                  to=c(2, 3, 5, 7, 4, 5, 7, 7, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17))

render_graph(create_graph(nodes_df=ndf) %>%
               add_edge_df(edf) %>%
               add_global_graph_attrs(attr='layout', value='dot', attr_type='graph') %>%
               add_global_graph_attrs(attr='fontcolor', value='black', attr_type='node') %>%
               add_global_graph_attrs(attr='shape', value='box', attr_type='node') %>%
               add_global_graph_attrs(attr='height', value=0.5, attr_type='node') %>%
               add_global_graph_attrs(attr='width', value=3, attr_type='node') %>%
               add_global_graph_attrs(attr='penwidth', value=2, attr_type='node') %>%
               add_global_graph_attrs(attr='fillcolor', value='white', attr_type='node') %>%
               select_nodes(conditions=type=='HMF') %>%
               set_node_attrs_ws(node_attr='color', value='#586BA4') %>%
               set_node_attrs_ws(node_attr='style', value='filled') %>%
               set_node_attrs_ws(node_attr='shape', value='diamond') %>%
               set_node_attrs_ws(node_attr='height', value=1) %>%
               clear_selection() %>%
               select_nodes(conditions=type=='TCGA') %>%
               set_node_attrs_ws(node_attr='shape', value='ellipse') %>%
               set_node_attrs_ws(node_attr='color', value='#F5DD90') %>%
               clear_selection() %>%
               select_nodes(conditions=type=='Breast') %>%
               set_node_attrs_ws(node_attr='shape', value='polygon') %>%
               set_node_attrs_ws(node_attr='sides', value=5) %>%
               set_node_attrs_ws(node_attr='height', value=1) %>%
               set_node_attrs_ws(node_attr='color', value='#688B58') %>%
               clear_selection() %>%
               select_nodes(conditions=type=='Ovarian') %>%
               set_node_attrs_ws(node_attr='shape', value='polygon') %>%
               set_node_attrs_ws(node_attr='sides', value=6) %>%
               set_node_attrs_ws(node_attr='color', value='#F68E5F')) %>%
  export_svg %>%
  charToRaw %>%
  rsvg_svg(file.path(supp_figs_dir, 'Phase_III_general_filtering.svg'))
