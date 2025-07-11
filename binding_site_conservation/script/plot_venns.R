# 加载必要的包
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(jsonlite)
library(VennDiagram)
library(optparse)
library(patchwork)
library(ggplot2)
library(grid)

# 定义命令行参数
option_list <- list(
  make_option(c("-j", "--json_file"), type = "character", default = "/data1/RiboGen/data/filter_specie_rbp.json",
              help = "JSON 文件路径，默认为 %default"),
  make_option(c("-o", "--output_file"), type = "character", default = "species_vs_human_venn.pdf",
              help = "Venn 图输出文件路径，默认为 %default")
)

opt <- parse_args(OptionParser(option_list = option_list))

# 检查 JSON 文件是否存在
if (file.exists(opt$json_file)) {
  # 读取 JSON 文件
  json_data <- read_json(opt$json_file)
  message("成功读取 JSON 文件")
  
  # 获取人类的 RBP 数据
  human_rbps <- json_data$human
  
  # 获取除人类外的其他物种名称
  other_species <- names(json_data)[names(json_data) != "human"]
  
  # 初始化存储 Venn 图的列表
  venn_plots <- list()
  
  # 遍历除人类外的每个物种
  for (i in seq_along(other_species)) {
    species <- other_species[i]
    species_rbps <- json_data[[species]]
    
    # 绘制 Venn 图
    venn_grob <- draw.pairwise.venn(
      area1 = length(human_rbps),
      area2 = length(species_rbps),
      cross.area = length(intersect(human_rbps, species_rbps)),
      category = c("Human", species),
      fill = c("lightblue", "lightgreen"),
      alpha = 0.5,
      cat.pos = c(0, 0),
      cat.dist = c(0.025, 0.025),
      scaled = FALSE,
      cat.fontface = 4,
      cat.fontfamily = "serif",
      cat.col = c("darkblue", "darkgreen"),
      fontfamily = "serif",
      lty = "blank",
      cex = 1.5
    )
    
    # 将 gList 转换为单个 grob 对象
    single_venn_grob <- grid::gTree(children = venn_grob)
    
    # 将 Venn 图转换为 ggplot 对象
    venn_ggplot <- ggplot() + 
      annotation_custom(grob = single_venn_grob) +
      theme_void() +
      theme(
        plot.margin = margin(10, 10, 10, 10),
        panel.spacing = unit(2, "cm")  # 增加子图之间的间距
      )
    
    venn_plots[[i]] <- venn_ggplot
  }
  
  # 动态计算列数，可根据实际情况调整
  ncol_value <- 2
  # 计算行数
  nrow_value <- ceiling(length(venn_plots) / ncol_value)
  
  # 使用 patchwork 拼图
  combined_plot <- wrap_plots(venn_plots, ncol = ncol_value, nrow = nrow_value) + 
    plot_layout(
      guides = "collect",
      widths = rep(1, ncol_value)  # 调整列宽
    ) +
    theme(
      panel.spacing.x = unit(3, "cm"),  # 增加列之间的间距
      panel.spacing.y = unit(1, "cm")   # 增加行之间的间距
    )
  
  # 设置 PDF 输出文件
  pdf(opt$output_file)
  
  # 输出拼图
  print(combined_plot)
  
  # 关闭 PDF 设备
  dev.off()
  message(sprintf("Venn 图已保存为 %s", opt$output_file))
} else {
  stop(sprintf("JSON 文件 %s 不存在", opt$json_file), call. = FALSE)
}
