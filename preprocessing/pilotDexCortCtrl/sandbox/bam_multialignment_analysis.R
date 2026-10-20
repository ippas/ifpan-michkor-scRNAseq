library(Rsamtools)
library(dplyr)
library(ggplot2)

bam_file <- "data/pilotDexCortCtrl_v2/BG60BLUE_L4-out/process/BG60BLUE_L4_merged_anno_barcode_headAligned.sorted.fixed.bam"
bam_file <- "data/pilotDexCortCtrl_v2/BG60BLUE_L4-out/process/BG60BLUE_L4_merged_anno_barcode_headAligned.sorted.bam"


# wczytaj interesujące tagi (NH, GN)
param <- ScanBamParam(tag=c("NH", "GN"))
param <- ScanBamParam(
  what = c("rname", "pos"),
  tag  = c("NH", "GN")
)

param <- ScanBamParam(
  what = c("rname", "pos", "mapq"),   # <--- tutaj dodaj mapq
  tag  = c("NH", "GN")
)

aln <- scanBam(bam_file, param=param)[[1]]

df <- data.frame(
  chr   = aln$rname,
  pos   = aln$pos,
  mapq  = aln$mapq,          # mapping quality
  gene  = unlist(aln$tag$GN),
  nh    = unlist(aln$tag$NH)
) %>%
  filter(gene != "")

df$mapq %>% count


aln <- scanBam(bam_file, param=param)[[1]]

nh <- aln$tag$NH
gn <- aln$tag$GN

df <- data.frame(
  gene = gn,
  nh = nh
) %>% 
  filter(gene!="") 

df %>% .$nh %>% count

df %>% filter(gene != "") %>% .$nh %>% unique()

# 1. Rozkład NH
nh_dist <- 
  df %>% .$nh %>% count

ggplot(nh_dist, aes(x=factor(x), y=freq)) +
  geom_col(fill="steelblue") +
  theme_minimal() +
  labs(title="Rozkład NH", x="NH", y="Liczba readów")

# 2. Multialignments per chromosome (na podstawie seqnames)
chr <- as.character(aln$rname)  # kolumna chromosomów
df_chr <- data.frame(chr=chr, nh=nh) %>%
  filter(nh > 1) %>%
  .$chr %>% count

ggplot(df_chr, aes(x=x, y=freq)) +
  geom_col(fill="darkred") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Multialignments per chromosome")

# 3. Multialignments per gene
df_gene <- df %>%
  filter(nh > 1) %>%head
  count(gene, name="multi_count")

# 4. Top 20 genów
top20 <- df_gene %>%
  arrange(desc(multi_count)) %>%
  head(20)

ggplot(top20, aes(x=reorder(gene, multi_count), y=multi_count)) +
  geom_col(fill="purple") +
  coord_flip() +
  theme_minimal() +
  labs(title="Top 20 genów z multialignments", x="Gen", y="Liczba readów")
