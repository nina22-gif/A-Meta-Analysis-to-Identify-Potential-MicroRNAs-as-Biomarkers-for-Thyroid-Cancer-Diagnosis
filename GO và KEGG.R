
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(patchwork)
library(dplyr)


#---------------------PTC-up-----------------------------------------------
mirna_list <- c("hsa-miR-146b-5p",
                "hsa-miR-221-3p",
                "hsa-miR-222-3p",
                "hsa-miR-21-5p")

all_targets <- data.frame()


for (mir in mirna_list) {
  result <- get_multimir(mirna = mir, table = "validated", summary = TRUE)
  all_targets <- rbind(all_targets, result@data)
}

# Lấy danh sách gen đích duy nhất
target_genes <- unique(all_targets$target_symbol)
write.csv(target_genes, "downregulated_miRNA_targets.csv", row.names = FALSE)

high_conf_targets <- all_targets %>%
  group_by(mature_mirna_id, target_symbol) %>%
  summarise(sources = n_distinct(database), .groups = "drop") %>%
  filter(sources >= 2)
# Lấy gene đích tin cậy
filtered_genes <- unique(high_conf_targets$target_symbol)
                         
#----------Phân tích GO và KEGG-------------
library(clusterProfiler)
library(org.Hs.eg.db)

# Chuyển SYMBOL → ENTREZID (để phân tích KEGG)
entrez_ids <- bitr(filtered_genes, fromType="SYMBOL",
                   toType="ENTREZID", OrgDb=org.Hs.eg.db)

# enrichGO cho BP, CC, MF
ego_bp <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05)

ego_cc <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05)

ego_mf <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05)

# enrichKEGG
ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                    organism = 'hsa', pvalueCutoff = 0.05)
p1 <- barplot(ego_bp, showCategory = 10, font.size = 12, title = "Biological Process (BP)")
p2 <- barplot(ego_cc, showCategory = 10, font.size = 12, title = "Cell Component (CC)")
p3 <- barplot(ego_mf, showCategory = 10, font.size = 12, title = "Molecular Function (MF)")
p4 <- barplot(ekegg, showCategory = 10, font.size = 12, title = "KEGG Pathway")


library(patchwork)
(p1 | p2) / (p3 | p4)

#------------PTC_down-----
# Danh sách miRNA của PTCdown
ptcdown_mirna <- c("hsa-miR-451a",
                   "hsa-miR-7-5p",
                   "hsa-miR-6749-5p",
                   "hsa-miR-6869-5p"
                  )

# Lấy tất cả gene đích đã được xác nhận (validated)
ptcdown_targets_all <- data.frame()
for (mir in ptcdown_mirna) {
  result_PTC_down <- get_multimir(mirna = mir, table = "validated", summary = TRUE)
  ptcdown_targets_all <- rbind(ptcdown_targets_all, result_PTC_down@data)
}

# Lưu danh sách tất cả gene đích
ptcdown_target_genes <- unique(ptcdown_targets_all$target_symbol)
write.csv(ptcdown_target_genes, "PTCdown_all_targets.csv", row.names = FALSE)


# Chuyển SYMBOL → ENTREZID
ptcdown_entrez <- bitr(ptcdown_target_genes, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Phân tích GO
ptcdown_go_bp <- enrichGO(gene = ptcdown_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05)

ptcdown_go_cc <- enrichGO(gene = ptcdown_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05)

ptcdown_go_mf <- enrichGO(gene = ptcdown_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.1)

# Phân tích KEGG
ptcdown_kegg <- enrichKEGG(gene = ptcdown_entrez$ENTREZID,
                           organism = 'hsa', pvalueCutoff = 0.05)

# Tạo biểu đồ
p1_ptcdown <- barplot(ptcdown_go_bp, showCategory = 10, font.size = 12, title = "Biological Process (PTCdown)")
p2_ptcdown <- barplot(ptcdown_go_cc, showCategory = 10, font.size = 12, title = "Cell Component (PTCdown)")
p3_ptcdown <- barplot(ptcdown_go_mf, showCategory = 10, font.size = 12, title = "Molecular Function (PTCdown)")
p4_ptcdown <- barplot(ptcdown_kegg, showCategory = 10, font.size = 12, title = "KEGG Pathway (PTCdown)")

# Ghép biểu đồ
(p1_ptcdown | p2_ptcdown) / (p3_ptcdown | p4_ptcdown)

#---------------MTC-down------------------------

# Danh sách miRNA MTC-down
mtcdown_mirna <- c("hsa-miR-451",
                   "hsa-miR-144-3p",
                   "hsa-miR-130a-3p",
                   "hsa-let-7i-5p"
                   )
# Lấy tất cả gene đích đã được xác nhận (validated)
mtcdown_targets_all <- data.frame()
for (mir in mtcdown_mirna) {
  result_MTC_down <- get_multimir(mirna = mir,table = "validated", summary = TRUE)
  mtcdown_targets_all <- rbind(mtcdown_targets_all,result_MTC_down@data)
}

# Lưu danh sách tất cả gene đích
mtcdown_target_genes <- unique(mtcdown_targets_all$target_symbol)
write.csv(mtcdown_target_genes, "MTCdown_all_targets.csv", row.names = FALSE)

# Chuyển SYMBOL → ENTREZID
mtcdown_entrez <- bitr(mtcdown_target_genes, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Phân tích GO
mtcdown_go_bp <- enrichGO(gene = mtcdown_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05)

mtcdown_go_cc <- enrichGO(gene = mtcdown_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05)

mtcdown_go_mf <- enrichGO(gene = mtcdown_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.1)

# Phân tích KEGG
mtcdown_kegg <- enrichKEGG(gene = mtcdown_entrez$ENTREZID,
                           organism = 'hsa', pvalueCutoff = 0.05)

# Tạo biểu đồ
p1_mtcdown <- barplot(mtcdown_go_bp, showCategory = 10, font.size = 12, title = "Biological Process (MTCdown)")
p2_mtcdown <- barplot(mtcdown_go_cc, showCategory = 10, font.size = 12, title = "Cell Component (MTCdown)")
p3_mtcdown <- barplot(mtcdown_go_mf, showCategory = 10, font.size = 12, title = "Molecular Function (MTCdown)")
p4_mtcdown <- barplot(mtcdown_kegg, showCategory = 10, font.size = 12, title = "KEGG Pathway (MTCdown)")

# Ghép biểu đồ
(p1_mtcdown | p2_mtcdown) / (p3_mtcdown | p4_mtcdown)

#------------------------MTC-UP--------------------------------
#------------------------MTC-UP--------------------------------
mtcup_mirna <- c("hsa-miR-7-5p",
                 "hsa-miR-129-3p",
                 "hsa-miR-376c-3p",
                 "hsa-miR-136-5p"
)

mtcup_targets_all <- data.frame()
for (mir in mtcup_mirna) {
  result_MTC_up <- get_multimir(mirna = mir, table = "validated", summary = TRUE)
  mtcup_targets_all <- rbind(mtcup_targets_all, result_MTC_up@data)
}

# Lưu danh sách tất cả gene đích
mtcup_target_genes <- unique(mtcup_targets_all$target_symbol)
write.csv(mtcup_target_genes, "MTCup_all_targets.csv", row.names = FALSE)

# Chuyển SYMBOL → ENTREZID
mtcup_entrez <- bitr(mtcup_target_genes, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Phân tích GO
mtcup_go_bp <- enrichGO(gene = mtcup_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05)

mtcup_go_cc <- enrichGO(gene = mtcup_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05)

mtcup_go_mf <- enrichGO(gene = mtcup_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.1)

# Phân tích KEGG
mtcup_kegg <- enrichKEGG(gene = mtcup_entrez$ENTREZID,
                         organism = 'hsa', pvalueCutoff = 0.05)

# Tạo biểu đồ
p1_mtcup <- barplot(mtcup_go_bp, showCategory = 10, font.size = 12, title = "Biological Process (MTCup)")
p2_mtcup <- barplot(mtcup_go_cc, showCategory = 10, font.size = 12, title = "Cell Component (MTCup)")
p3_mtcup <- barplot(mtcup_go_mf, showCategory = 10, font.size = 12, title = "Molecular Function (MTCup)")
p4_mtcup <- barplot(mtcup_kegg, showCategory = 10, font.size = 12, title = "KEGG Pathway (MTCup)")

# Ghép biểu đồ
(p1_mtcup | p2_mtcup) / (p3_mtcup | p4_mtcup)
#---------------------Chuan bi du lieu ve CYTOSCAPE------------------
#----------------PTC---------------------

all_target_genes_PTC <- unique(c(filtered_genes, ptcdown_target_genes))  # hoặc ptcdown_filtered_genes nếu bạn lọc theo nguồn

# Chuyển SYMBOL → ENTREZID
entrez_df_PTC <- bitr(all_target_genes_PTC, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 2. Tạo đối tượng STRINGdb
string_db_PTC <- STRINGdb$new(version = "11.5", species = 9606,
                              score_threshold = 700, input_directory = "")

# 3. Ánh xạ ENTREZID → STRING ID
mapped_genes_PTC <- string_db_PTC$map(entrez_df_PTC, "ENTREZID", removeUnmappedRows = TRUE)

# 4. Lấy dữ liệu tương tác
interactions_PTC <- string_db_PTC$get_interactions(mapped_genes_PTC$STRING_id)

# 5. Lọc tương tác có độ tin cậy cao (score ≥ 950)
filtered_interactions_PTC <- interactions_PTC %>%
  filter(combined_score >= 950)

nrow(filtered_interactions_PTC)  # Kiểm tra số lượng

# 6. Ánh xạ ID → tên gen (SYMBOL)
alias_map_PTC <- mapped_genes_PTC[, c("STRING_id", "SYMBOL")]

# Đổi ID 'from' thành tên gen
inter1_PTC <- merge(filtered_interactions_PTC, alias_map_PTC,
                    by.x = "from", by.y = "STRING_id", all.x = TRUE)
colnames(inter1_PTC)[colnames(inter1_PTC) == "SYMBOL"] <- "from_name"

# Đổi ID 'to' thành tên gen
inter2_PTC <- merge(inter1_PTC, alias_map_PTC,
                    by.x = "to", by.y = "STRING_id", all.x = TRUE)
colnames(inter2_PTC)[colnames(inter2_PTC) == "SYMBOL"] <- "to_name"

# 7. Chọn cột cần thiết và đổi tên
ppi_named_PTC <- inter2_PTC[, c("from_name", "to_name", "combined_score")]
colnames(ppi_named_PTC) <- c("Source", "Target", "Score")

# 8. Xuất file cho Cytoscape
write.table(ppi_named_PTC,
            file = "PPI_network_PTCmiRNA.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#------------------MTC-----------------
unique_targets_MTCup <- mtcup_targets_all %>%
  distinct(mature_mirna_id, target_symbol, .keep_all = TRUE)

unique_targets_MTCdowwn <- mtcdown_targets_all %>%
  distinct(mature_mirna_id, target_symbol, .keep_all = TRUE)

all_target_genes_MTC <- unique(c(unique_targets_MTCup$target_symbol, unique_targets_MTCdowwn$target_symbol))

# Chuyển SYMBOL → ENTREZID
entrez_df_MTC <- bitr(all_target_genes_MTC, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 2. Tạo đối tượng STRINGdb
string_db_MTC <- STRINGdb$new(version = "11.5", species = 9606,
                              score_threshold = 700, input_directory = "")

# 3. Ánh xạ ENTREZID → STRING ID
mapped_genes_MTC <- string_db_MTC$map(entrez_df_MTC, "ENTREZID", removeUnmappedRows = TRUE)

# 4. Lấy dữ liệu tương tác
interactions_MTC <- string_db_MTC$get_interactions(mapped_genes_MTC$STRING_id)

# 5. Lọc tương tác có độ tin cậy cao (score ≥ 950)
filtered_interactions_MTC <- interactions_MTC %>%
  filter(combined_score >= 950)

nrow(filtered_interactions_MTC) # Kiểm tra số lượng

# 6. Ánh xạ ID → tên gen (SYMBOL)
alias_map_MTC <- mapped_genes_MTC[, c("STRING_id", "SYMBOL")]

# Đổi ID 'from' thành tên gen
inter1_MTC <- merge(filtered_interactions_MTC, alias_map_MTC,
                    by.x = "from", by.y = "STRING_id", all.x = TRUE)
colnames(inter1_MTC)[colnames(inter1_MTC) == "SYMBOL"] <- "from_name"

# Đổi ID 'to' thành tên gen
inter2_MTC <- merge(inter1_MTC, alias_map_MTC,
                    by.x = "to", by.y = "STRING_id", all.x = TRUE)
colnames(inter2_MTC)[colnames(inter2_MTC) == "SYMBOL"] <- "to_name"

# 7. Chọn cột cần thiết và đổi tên
ppi_named_MTC <- inter2_MTC[, c("from_name", "to_name", "combined_score")]
colnames(ppi_named_MTC) <- c("Source", "Target", "Score")

# Chuẩn hóa chiều tương tác (A-B = B-A)
ppi_named_MTC_unique <- ppi_named_MTC %>%
  rowwise() %>%
  mutate(
    gene1 = min(Source, Target),
    gene2 = max(Source, Target)
  ) %>%
  ungroup() %>%
  distinct(gene1, gene2, .keep_all = TRUE) %>%
  select(Source = gene1, Target = gene2, Score)

# 8. Xuất file cho Cytoscape
write.table(ppi_named_MTC_unique,
            file = "PPI_network_MTCmiRNA_unique.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)















