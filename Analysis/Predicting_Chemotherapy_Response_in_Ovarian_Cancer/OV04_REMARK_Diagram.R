# Clean environment
freshr::freshr()

# Libraries
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(dplyr)
library(this.path)

# Paths
pancancer_dir <- dirname(this.path())
base_dir <- dirname(dirname(pancancer_dir))
figs_dir <- file.path(base_dir, 'Figures')
supp_figs_dir <- file.path(figs_dir, 'Supplementary')
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_figs_dir, showWarnings = FALSE, recursive = TRUE)

# Render graph
ndf_a <- data.frame(id=1:11, type=c('Original', rep('Filters', 7), rep('Filtered', 3)),
                  label=c("Primary Patient\nCohort\n\n130 Patients",
                          "Remove Patients\nWith Low Quality\nSequencing Data\n\n33 Patients",
                          "Remove Patients\nWithout Treatment\nData\n\n11 Patients",
                          "Remove Patients\nWithout CA125\nData\n\n17 Patients",
                          "Remove Patients\nWith Gaps in\nCA125 Data\n\n18 Patients",
                          "Remove Patients\nWithout Abnormal\n1st-line CA125\n\n2 Patients",
                          "Remove Patients\nWith Missing\nCovariates\n\n2 Patients",
                          "Remove Patients\nEnrolled in Clinical Trials\n\n2 Patients",
                          "Curated Cohort for\nPlatinum Analysis\n\n45 Patients", 
                          "Matched Tissue-\nTSO500 Cohort\n\n8 Patients",
                          "Matched Tissue-\nPlasma Cohort\n\n9 Patients"))
edf_a <- data.frame(id=1:10,
                  from=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9),
                  to=c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

ndf_b <- data.frame(id=1:7, type=c('Filtered', 'Filters', 'Original', 'Filters', rep('Filtered', 3)),
                    label=c("Curated Cohort for\nPlatinum Analysis\n\n45 Patients",
                            "Remove Patients\nNot Receiving\nTaxanes\n\n21 Patients",
                            "Primary Taxane\nCohort\n\n24 Patients",
                            "Add Patients With\nCA125 Gaps\nin Platinum\n\n5 Patients",
                            "Curated Cohort for\nTaxane Analysis\n\n29 Patients",
                            "Matched Tissue-\nTSO500 Cohort\n\n3 Patients",
                            "Matched Tissue-\nPlasma Cohort\n\n4 Patients"))
edf_b <- data.frame(id=1:6,
                    from=c(1, 2, 3, 4, 5, 5),
                    to=c(2, 3, 4, 5, 6, 7))

ndf_c <- data.frame(id=1:12, type=c('Original', rep('Filters', 2), 'Filtered', 'Original', 'Filters', rep('Filtered', 2), 'Filters', rep('Filtered', 3)),
                    label=c("Organoid Samples\n\n11 Patients",
                            "Remove Noisy\nSamples\n\n1 Patient",
                            "Remove Patients\nWith Spheroids\n\n2 Patients",
                            "Final Organoid\nCohort\n\n8 Patients",
                            "Spheroid Samples\n\n19 Patients",
                            "Remove Noisy\nSamples\n\n4 Patients",
                            "Final Spheroid\nCohort\n\n15 Patients",
                            "Curated Cohort for\nPlatinum Analysis\n\n45 Patients",
                            "Remove Patients\nNot Receiving\nAnthracyclines\n\n15 Patients",
                            "Curated Cohort for\nAnthracycline\nAnalysis\n\n30 Patients",
                            "Matched Tissue-\nTSO500 Cohort\n\n2 Patients",
                            "Matched Tissue-\nPlasma Cohort\n\n5 Patients"))
edf_c <- data.frame(id=1:9,
                    from=c(1, 2, 3, 5, 6, 8, 9, 10, 10),
                    to=c(2, 3, 4, 6, 7, 9, 10, 11, 12))


render_graph(create_graph(nodes_df=ndf_a) %>%
               add_edge_df(edf_a) %>%
               add_global_graph_attrs(attr='layout', value='dot', attr_type='graph') %>%
               add_global_graph_attrs(attr='fontcolor', value='black', attr_type='node') %>%
               add_global_graph_attrs(attr='shape', value='box', attr_type='node') %>%
               add_global_graph_attrs(attr='height', value=0.8, attr_type='node') %>%
               add_global_graph_attrs(attr='width', value=1.75, attr_type='node') %>%
               add_global_graph_attrs(attr='penwidth', value=2, attr_type='node') %>%
               add_global_graph_attrs(attr='fillcolor', value='white', attr_type='node') %>%
               select_nodes(conditions=type=='Filters') %>%
               set_node_attrs_ws(node_attr='fillcolor', value='#E6E6E6') %>%
               clear_selection() %>%
               select_nodes(conditions=type=='Filtered') %>%
               set_node_attrs_ws(node_attr='fillcolor', value='#D9FFFF') %>%
               clear_selection()) %>%
    export_svg %>%
    charToRaw %>%
    rsvg_svg(file.path(supp_figs_dir, 'OV04_REMARK_Platinum.svg'))

render_graph(create_graph(nodes_df=ndf_b) %>%
               add_edge_df(edf_b) %>%
               add_global_graph_attrs(attr='layout', value='dot', attr_type='graph') %>%
               add_global_graph_attrs(attr='fontcolor', value='black', attr_type='node') %>%
               add_global_graph_attrs(attr='shape', value='box', attr_type='node') %>%
               add_global_graph_attrs(attr='height', value=0.8, attr_type='node') %>%
               add_global_graph_attrs(attr='width', value=1.75, attr_type='node') %>%
               add_global_graph_attrs(attr='penwidth', value=2, attr_type='node') %>%
               add_global_graph_attrs(attr='fillcolor', value='white', attr_type='node') %>%
               select_nodes(conditions=type=='Filters') %>%
               set_node_attrs_ws(node_attr='fillcolor', value='#E6E6E6') %>%
               clear_selection() %>%
               select_nodes(conditions=type=='Filtered') %>%
               set_node_attrs_ws(node_attr='fillcolor', value='#D9FFFF') %>%
               clear_selection()) %>%
    export_svg %>%
    charToRaw %>%
    rsvg_svg(file.path(supp_figs_dir, 'OV04_REMARK_Taxane.svg'))

render_graph(create_graph(nodes_df=ndf_c) %>%
               add_edge_df(edf_c) %>%
               add_global_graph_attrs(attr='layout', value='dot', attr_type='graph') %>%
               add_global_graph_attrs(attr='fontcolor', value='black', attr_type='node') %>%
               add_global_graph_attrs(attr='shape', value='box', attr_type='node') %>%
               add_global_graph_attrs(attr='height', value=0.8, attr_type='node') %>%
               add_global_graph_attrs(attr='width', value=1.75, attr_type='node') %>%
               add_global_graph_attrs(attr='penwidth', value=2, attr_type='node') %>%
               add_global_graph_attrs(attr='fillcolor', value='white', attr_type='node') %>%
               select_nodes(conditions=type=='Filters') %>%
               set_node_attrs_ws(node_attr='fillcolor', value='#E6E6E6') %>%
               clear_selection() %>%
               select_nodes(conditions=type=='Filtered') %>%
               set_node_attrs_ws(node_attr='fillcolor', value='#D9FFFF') %>%
               clear_selection()) %>%
    export_svg %>%
    charToRaw %>%
    rsvg_svg(file.path(supp_figs_dir, 'OV04_REMARK_Anthracycline.svg'))
