#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tools")
library("tidyr")
library("xlsx")

# Parameter validation
args_parser = ArgumentParser(
  prog="selex_rpm_top_n",
  description=""
)
args_parser$add_argument("--in-rpm-csv", "-i", help="Input file must be csv as created by SELEXrpm", required=TRUE)
args_parser$add_argument("--out-file", "-o", default="top1000.xlsx", help="Output file. Default: top1000.xlsx")
args_parser$add_argument("--top", "-n", type="integer", default=1000, help="Output top n sequences from every round.")
#args <- args_parser$parse_args(c("-i", "rpm.csv"))
args <- args_parser$parse_args()

# EXCEL Workbook
wb <- createWorkbook(type="xlsx")

style_title = CellStyle(wb) + Font(wb, isBold=TRUE, heightInPoints = 11)
style_text = CellStyle(wb) + Font(wb, isBold=TRUE, heightInPoints = 10)
style_colnames = CellStyle(wb) + Font(wb, isBold=TRUE, heightInPoints = 10)
style_data = CellStyle(wb) + Font(wb, name="Courier New", heightInPoints = 10)
style_number = CellStyle(wb, dataFormat=DataFormat("0")) + Font(wb, heightInPoints = 10) # no decimal point

# Data input
df_rpm <-  read.csv(args$in_rpm_csv, sep=";", header=TRUE)
round_names <- levels(factor(df_rpm$round))
for (round_name in round_names) {
  # Data stuff
  df_round <- df_rpm %>%
    filter(round == round_name) %>%
    arrange(desc(rpm)) %>%
    mutate(rank = row_number(desc(rpm))) %>% # rpm=round(rpm, digits=0)
    select(rank, count, rpm, seq, p5_seq_p3)
  
  df_top_1000 <- df_round %>%
    top_n(-args$top, rank) %>% # select from bottom by rank
    rename("absolute"="count", "rpm"="rpm", "random region"="seq", "full sequence"="p5_seq_p3")
  df_top_1000$"analysed as" <- " "
  
  # Fill sheet with info
  sheet <- createSheet(wb, sheetName=round_name)
  title_row <- createRow(sheet, rowIndex=1)
  title_cell <- createCell(title_row, colIndex=1)
  setCellValue(title_cell[[1,1]], paste("Top ", args$top, " reads, Illumina, Round ", round_name))
  setCellStyle(title_cell[[1,1]], style_title)
  
  empty_row = createRow(sheet, rowIndex=2)
  
  round_specific_row = createRow(sheet, rowIndex = 3)
  round_name_cell = createCell(round_specific_row, colIndex = 1)
  round_total_cell = createCell(round_specific_row, colIndex = 2)
  
  setCellValue(round_name_cell[[1,1]], round_name)
  setCellStyle(round_name_cell[[1,1]], style_text)
  setCellValue(round_total_cell[[1,1]], paste(sum(df_round$count), "reads", sep=" "))
  setCellStyle(round_total_cell[[1,1]], style_text)
  
  setColumnWidth(sheet, 1, 7)
  setColumnWidth(sheet, 2, 8)
  setColumnWidth(sheet, 3, 8)
  setColumnWidth(sheet, 4, 40)
  setColumnWidth(sheet, 5, 15)
  
  addDataFrame(
    df_top_1000, 
    sheet,
    col.names = TRUE,
    row.names = FALSE,
    colnamesStyle = style_colnames,
    startRow=4,
    startCol=1,
    colStyle = list(`3`=style_number, `4`=style_data, `5`=style_data)
  )

}

saveWorkbook(wb, args$out_file)
