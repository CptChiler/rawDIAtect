
# set path to DIA-NN quant folder (multiple quants possbile, only .tsv in this folder, subfolders are allowed)
path_in = "C:/Users/blumenscheitc/Desktop/quants_sorting/full"

# filter (0 = off)
pep_filter = 3


# whiteliste (1 = on)
whitelist = 1

# Set Experiment Name (if not leave blank "", no space!) 
Exp_name = "only_PCR_hits"

# Set paths for Report (Will create an Folder /Report_*Experiment_Name*)
# if not leave blank "" , no space! Output will be generated in input folder!
path_out = ""

# Set path to aro index database (*.tsv)
path_db_in = "../database/aro_index_DIA.tsv"


# check Input and Output this scipt overwrites old results!
# Close opened pdf Results or it will crash.
# run Source or Ctrl+Shift+S 




# Dont change things below this !
#dependencies#############################################################################################################

options(warn=-1)
if(!is.null(dev.list())) dev.off()
suppressMessages(
  
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos='https://cloud.r-project.org/')
  sapply(pkg, require, character.only = TRUE)
}
)

# Usage example
packages<-c("readr","tidyverse","profvis","tidyr", "ggplot2", "dplyr", "cowplot", "plyr", "devtools", "gameofthrones","tm","crayon","progress","Hmisc")
suppressMessages(
  check.packages(packages)
  )
print("Checking if dependencies are available =")
print(check.packages(packages))


#rawSplit################################################################################################################

cat("\n")
cat(bold(underline("Loading module rawSPLIT")))
cat("\n")
cat("\n")

# gets filenames_tsv from path_in
filenames_tsv <- list.files(path=path_in,
                        pattern=".tsv")

l1 = length(filenames_tsv)
l2 = paste("Found",l1,sep =" ")
l3 = paste(l2,"DIA-NN Quants in")
cat(blue(paste(l3, path_in,sep=" ")) %+%
      "\n") 
cat("\n")

for(i in filenames_tsv) 
  {
  filepath <- file.path(path_in,paste(i,sep=""))
  t1 = paste("Loading DIA-NN", i,sep = " ")
  t2 = paste(t1, "from",sep = " ")
  print(paste(t2,path_in,sep = " "))
  temp_dat <- read_delim(filepath, "\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols())
  cat(bgGreen$black$bold(paste("Starts splitting", i, sep= " ")) %+%
        "\n")
  
  #read quant
  quant <- temp_dat
  
  
  #Spliting
  Sample_raw <- basename(as.character(quant$File.Name))
  Sample <- removeWords(Sample_raw, ".raw")
  Samples <- as.data.frame(Sample)
  Samples$rownumber = 1:nrow(Samples)
  quant$rownumber = 1:nrow(quant)
  quant_split <- merge(quant,Samples, by="rownumber")
  quant_split <- quant_split %>% select(-rownumber)
  
  # number of Files
  
  file_count <- length(unique(quant_split$Sample))
  
  file_print <- paste(file_count,"splitted Quant files saved to", sep = " ")
  
  #Create Dir
  
   out_dir <- paste(path_in,"split/", sep = "/")
      dir.create(path_out, showWarnings = FALSE)
      dir.create(out_dir, showWarnings = FALSE)
  
  # Split by variable
  spt2 <- split(quant_split, quant_split$Sample)
  
  
  pb <- progress_bar$new(format = " progress [:bar] :percent eta: :eta", # add
                         total = length(names(spt2)), clear = FALSE, width= 60)
  
  # check if output file are present
  list_new <- names(spt2)
  list_new_csv=c(paste0(list_new, ".csv"))
  list_old_csv <- list.files(path=out_dir,
                              pattern=".csv")
  matches <- list_new_csv %in% list_old_csv
  df_matches <- data.frame(list_new,matches)
  names(df_matches)[1] <- "names"
  df_matches_filt <- filter(df_matches, matches == "FALSE")
    list_needed_csv <- as.character(df_matches_filt$names)

  # Save
  lapply(list_needed_csv, function(x){
    pb$tick()
    write_csv(spt2[[x]], path = paste(out_dir,x, "_quant.csv", sep = ""))
    })
  
  # loop end msg
  cat(underline(paste("Done splitting with", i, sep=" ") %+%
                  "\n" %+%
                  paste(file_print,out_dir,sep = " "))%+%
        "\n")
  cat("\n")
  cat("\n")

}

cat(paste("Saved all splitted Quant files in", out_dir,sep = " ")%+%
      "\n" )

#rawCRIT############################################################################

cat("\n")
cat(bold(underline("Loading module CRIT")))
cat("\n")
cat("\n")

# gets filenames from path_in
path_split <- paste(path_in,"split",sep = "/")
filenames <- list.files(path=path_split,
                        pattern=".csv")

cat(green(
  paste("Found splitted Quants in", path_in,sep = " ") %+%
    "\n" %+%
    paste("Found Index in", path_db_in,sep = " "))%+%
    "\n")

# set to clear theme
theme_set(theme_cowplot())

# reads ARO-Index
aro_index_DIA <- read_delim(path_db_in, "\t", escape_double = FALSE, trim_ws = TRUE,col_types = cols())

cat(underline(green("ARO Index was successfully loaded!")) %+%
      "\n" %+%
      "\n" %+%
      "\n")

# starts popup progress
pb <- winProgressBar(title="Total Progress of rawDetect", label="0% done", min=0, max=100, initial=0)
loop_count <- 0


cat(bgCyan(bold("Starting analyzing with following options:","\n")))
cat(bold(paste("Unique Peptide filter = ", pep_filter)),"\n","\n","\n")

# main loops start
for(i in filenames){
  
  # get file path from file.names path.in
  filepath <- file.path(path_split,paste(i,sep=""))
  
  cat(bold(paste("Loading",i,sep = " "))%+%
        "\n") 
  
  # input from loop (filepath)
  temp_dat <-  read_delim(filepath, ",", escape_double = FALSE, trim_ws = TRUE,col_types = cols())
  
  # paste start
  
  cat(bgGreen$black$bold(head(paste("Starts analyzing",head(paste(temp_dat$Sample),1),sep = " "),1)) %+%
        "\n" %+%
        "\n")
  
  # greps only ARO tags from temp
  quant_ARO <- temp_dat[grep("ARO:", temp_dat$Genes), ]
  
  # only selects (Genes, Stripped.Sequence, Precursor.Quantity)
  quant_ARO_small <- select(quant_ARO,c("Sample","Genes","Stripped.Sequence","Precursor.Quantity"))
  
  # sorts PQ.seq and PQ.quant (highest top)
  quant_ARO_small_sort <- quant_ARO_small[order(quant_ARO_small$Stripped.Sequence,-abs(quant_ARO_small$Precursor.Quantity)),]
  
  # removes duplicate PQ.seq keeps highest PQ.quant
  quant_ARO_small_sort_uni <- quant_ARO_small_sort[ !duplicated(quant_ARO_small_sort$Stripped.Sequence), ] 
  # input for count amr
  
  # seperate multiple aro hits
  quant_ARO_small_sort_uni_sep <- separate_rows(quant_ARO_small_sort_uni,Genes,sep=";",convert = TRUE)
  quant_ARO_sep <- quant_ARO_small_sort_uni_sep[grep("ARO:", quant_ARO_small_sort_uni_sep$Genes), ]
  # input for isoformes
  
  # count amr_gene_family
  quant_ARO_sep_uni_pep <- quant_ARO_sep[ !duplicated(quant_ARO_sep$Stripped.Sequence), ] 
  quant_ARO_AMR_index<- merge(x = quant_ARO_sep_uni_pep, y = aro_index_DIA[ , c("Genes", "amr_gene_family")], by = "Genes", all.x=TRUE)
  quant_ARO_count <- quant_ARO_AMR_index %>% group_by(amr_gene_family) %>% add_count(amr_gene_family,name = "peptides")
  
  if(nrow(quant_ARO_count) == 0) {
  mock_data <- data.frame(matrix(ncol = 9, nrow = 0))
  col_name_mock <- c("Sample","Names","Precursor.Quantity.Top3","amr_gene_family","Subfamily","peptides","resistance_mechanism","drug_class","Note")
  colnames(mock_data) <- col_name_mock
  mock_data[nrow(mock_data) + 1,] = c("None",head(paste(temp_dat$Sample),1),"None",0,0,"None","None","None","None")
  ggplot_filter_small_drugs <- mock_data
  } else {
  # top3 calculation (1.group/2.take top3/3.mean group/4.merge by amr_group)
  quant_top3 <- quant_ARO_count %>%
    group_by(amr_gene_family) %>%
    top_n(n = 3, wt = Precursor.Quantity) #Top3
  quant_top3_mean <- aggregate(quant_top3[, 4], list(quant_top3$amr_gene_family), mean)  #mean of group
  names(quant_top3_mean)[1] <- "amr_gene_family"
  quant_ARO_index_top3 <- merge(x = quant_ARO_count, y = quant_top3_mean[ , c("amr_gene_family", "Precursor.Quantity")], by = "amr_gene_family", all.x=TRUE)
  
  # sort data 
  quant_ARO_index_top3_sort <- select(quant_ARO_index_top3,c("Sample","Stripped.Sequence","Precursor.Quantity.y","amr_gene_family","peptides"))
  names(quant_ARO_index_top3_sort)[3] <- "Precursor.Quantity.Top3"
  
  if (pep_filter > 1) 
  {
    quant_ARO_index_top3_filt <- quant_ARO_index_top3_sort %>% filter(peptides > pep_filter)
  }
  
  # gets iso_names
  quant_ARO_iso_index <- merge(x = quant_ARO_sep, y = aro_index_DIA[ , c("Genes","Names","amr_gene_family","Subfamily","Note")], by = "Genes", all.x=TRUE)
  quant_ARO_iso_count <- quant_ARO_iso_index %>% group_by(amr_gene_family) %>% add_count(Names,name = "count_name")
  
  max_isos <- summarize(quant_ARO_iso_count$count_name,quant_ARO_iso_count$amr_gene_family,max)
    names(max_isos)[1] <- "amr_gene_family"
    names(max_isos)[2] <- "max_count_name"
  quant_ARO_iso_count_max <- merge(quant_ARO_iso_count,max_isos)
  quant_ARO_iso_count_diff <- quant_ARO_iso_count_max %>% mutate(diff = max_count_name-count_name)
  iso_name_best <- quant_ARO_iso_count_diff %>% filter(diff <= 2)
  iso_name_best_small <- select(iso_name_best,c("amr_gene_family","Genes","Names","Subfamily","Note"))
  iso_name_best_small_uni <-  unique(iso_name_best_small)
  
  # output with corrected isoform
  quant_ARO_name <- merge(quant_ARO_index_top3_filt,iso_name_best_small_uni)
  quant_ARO_iso <- left_join(x = quant_ARO_name, y = aro_index_DIA[ , c("Genes","resistance_mechanism","drug_class")], by = "Genes")
  
  # whitelisting for resistance_mechanism and split for ggplot
  if(whitelist == 1) {patterns <- c("kpc","oxa","shv","tem","van","mcr","cmy","aac","ctx","ndm","vim")
  ggplot_filter <- filter(quant_ARO_iso, grepl(paste(patterns, collapse="|"), Names,ignore.case = TRUE))
  ggplot_filter_small <- unique(select(ggplot_filter,c("Sample","Names","Precursor.Quantity.Top3","amr_gene_family","Subfamily","peptides","resistance_mechanism","drug_class","Note")))
  ggplot_filter_small_drugs <- unique(separate_rows(ggplot_filter_small,drug_class,sep=";",convert = TRUE))} else 
  {ggplot_filter_small_drugs <- quant_ARO_iso}
  
  }
  
  # if no pattern found print mock found for visulication
  if(nrow(ggplot_filter_small_drugs) == 0) {
    ggplot_filter_small_drugs[nrow(ggplot_filter_small_drugs) + 1,] = c("None",head(paste(temp_dat$Sample),1),"None",0,0,"None","None","None","None","None","None")
    ggplot_input <- na.omit(ggplot_filter_small_drugs)
  } else {
    ggplot_input <- ggplot_filter_small_drugs
  }
  
  #data
  #ggplot_input = visulize
  #quant_ARO_iso = table_report
  
  # visulize data
  
  # predefine colors
  
 if(length(unique(ggplot_input$amr_gene_family)) > 13) 
   {group.fills <- got(length(unique(ggplot_input$amr_gene_family)), option = "Daenerys")}
  else {  group.fills <-  c("AAC(6')" = "#541a1c",
                            "KPC beta-lactamase" = "#6C3237FF",
                            "OXA beta-lactamase" = "#604147FF",
                            "SHV beta-lactamase" = "#545058FF",
                            "TEM beta-lactamase" = "#46606AFF",
                            "CMY beta-lactamase" = "#38707CFF",
                            "CTX-M beta-lactamase" = "#2B818EFF",
                            "NDM beta-lactamase" = "#478B91FF",
                            "VIM beta-lactamase" = "#639594FF",
                            "MCR phosphoethanolamine transferase" = "#80A098FF",
                            "van ligase" = "#9AA99BFF",
                            "vanR" = "#B5B39EFF",
                            "vanS" = "#D1BDA2FF")
  }

  group.colors <- c(broadspectrum = "#006400",
                    Carbapenemase = "#DC143C",
                    ESBL ="#ffcccb",
                    Narrow = "#90ee90",
                    narrowspectrum = "#32CD32")

  # unique peptide count
  
  ggplot_pep_count <- unique(select(ggplot_input,c("Sample","amr_gene_family","peptides")))
  ggplot_pep_names <- unique(select(ggplot_input,c("Sample","amr_gene_family","Names","Subfamily","Note")))
  ggplot_pep_count_names <- ggplot_pep_names %>% group_by(Names) %>% add_count(Names,name = "count_1")
  count_amr_name <- count(ggplot_pep_count_names$amr_gene_family)
  names(count_amr_name)[1] <- "amr_gene_family"
  ggplot_pep_count2_names <- merge(ggplot_pep_count_names,count_amr_name)
  ggplot_names <- ggplot_pep_count2_names %>% mutate(perc = round((count_1/freq)*100,1))  
  
  p1=
    ggplot() +
    geom_bar(data = ggplot_pep_count,
             aes(
               x = reorder(amr_gene_family,-as.numeric(peptides)),
               y = as.numeric(peptides),
               fill = amr_gene_family
             ),
             stat="identity") +
    labs(x = "AMR Family",
         y = "Count of unique Hits") +
    geom_text(data = ggplot_pep_count,
              aes(
                x = reorder(amr_gene_family,-as.numeric(peptides)),
                y = as.numeric(peptides),
                label = peptides),
              vjust = 2) +
    facet_grid(~Sample) +
    theme(plot.subtitle = element_text(vjust = 1), 
          plot.caption = element_text(vjust = 1), 
          axis.text = element_text(hjust = 0), 
          axis.text.x = element_text(angle = -55,
                                     size = 9),
          legend.position = "none") +
    scale_fill_manual(values = group.fills,
                      na.value="white",na.translate=FALSE)
  
  pie_2=
    ggplot(data = ggplot_names) +
    geom_bar(aes(
      x = Sample,
      y = perc,
      fill = amr_gene_family
    ),
    stat="identity") + 
    geom_text(aes(
      x = Sample,
      y = perc,
      label = Names
    ),
    position = position_stack(vjust = 1))+
    geom_label(aes(x = Sample,
                   y = perc,
                   label = Subfamily),
               position = position_stack(vjust = 1))+
    facet_wrap(~amr_gene_family) +
    scale_fill_manual(values = group.fills,
                      na.value="white",na.translate=FALSE) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none")
  
  sup_pie_2 <- pie_2 + coord_polar("y", start=0)
  

  
  # Top3 per AMR
  
  ggplot_PQ_amr <- unique(select(ggplot_input,c("Sample","amr_gene_family","Precursor.Quantity.Top3")))
  
  p2=
    ggplot(data = ggplot_PQ_amr) +
    geom_bar(
      aes(
        x = reorder(amr_gene_family,-as.numeric(Precursor.Quantity.Top3)),
        y = as.numeric(Precursor.Quantity.Top3),
        fill = amr_gene_family
      ),
      stat="identity"
    ) +
    labs(x = "AMR Family",
         y = "Precursor Quantity (Top3)") +
    facet_grid(~Sample) +
    theme(plot.subtitle = element_text(vjust = 1), 
          plot.caption = element_text(vjust = 1), 
          axis.text = element_text(hjust = 0), 
          axis.text.x = element_text(angle = -55,
                                     size = 9),
          legend.position = "none") +
    scale_fill_manual(values = group.fills,
                      na.value="white",na.translate=FALSE)
  
  # gets matrix for main_page
  amr_gene_family <- c("AAC(6')","KPC beta-lactamase","OXA beta-lactamase","SHV beta-lactamase","TEM beta-lactamase","CMY beta-lactamase","CTX-M beta-lactamase",
                       "NDM beta-lactamase","VIM beta-lactamase","MCR phosphoethanolamine transferase","van ligase","vanR","vanS")
  Found_AMR <- c("")
  main_page_dat <- data.frame(amr_gene_family,Found_AMR)
  
  # merging data for plot for main page
  suppressMessages(
    main_page_pep <- left_join(main_page_dat,ggplot_pep_count)
    )
  suppressMessages(
    main_page_pep_PQ <- left_join(main_page_pep,ggplot_PQ_amr)
    )
  suppressMessages(
    main_page_pep_PQ_names <- left_join(main_page_pep_PQ,ggplot_names)
    )

  main_page_pep_PQ_names <- main_page_pep_PQ_names %>% mutate(Found_AMR = ifelse(Sample =="","", amr_gene_family))

  # new main page
  
  main_page <- ggplot(data = main_page_pep_PQ_names) +
    geom_tile(aes(x = amr_gene_family,
                  y = Note,
                  fill = Found_AMR,
                  color = Note),
              size =2) +
    theme(plot.subtitle = element_text(vjust = 1), 
          plot.caption = element_text(vjust = 1), 
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          axis.text.x = element_text(angle = 45,
                                     size = 9,
                                     hjust = 1)) +
    geom_label(aes(x = amr_gene_family,
                  y = Note,
                  label = Names),
              vjust =-0.8) +
    geom_label(aes(x = amr_gene_family,
                   y = Note,
                   label = Subfamily),
               vjust =-2) +
    geom_label(aes(x = amr_gene_family,
                  y = Note,
                  label = peptides),
              vjust =1) +
    scale_fill_manual(values = group.fills,
                      na.value="white",na.translate=FALSE)  +
    scale_color_manual(values=group.colors,
                         na.value="white",na.translate=FALSE) +
    guides(fill = F) +
    facet_grid(~head(temp_dat$Sample,1))
  
  
  # drug_class and spectrum
  
  ggplot_drug_class <- select(ggplot_input,c("Sample","drug_class","peptides"))
  ggplot_drug_class <- unique(separate_rows(ggplot_drug_class,drug_class,sep=";",convert = TRUE))
  
  
  # 
  if(ggplot_drug_class$drug_class == "None" ) {
    
    page1 <- p1
    page2 <- p2
    
  } else {
    # gets sums and percent of drug_classes
    
    ggplot_drug_sum <- aggregate(ggplot_drug_class$peptides, by=list(drug_class=ggplot_drug_class$drug_class), FUN=sum)
    ggplot_drug_perc <- ggplot_drug_sum %>% mutate(perc = round((x/sum(x))*100,2))
    ggplot_drug_perc$Sample <- rep(head(ggplot_drug_class$Sample,1),nrow(ggplot_drug_perc)) # make new column 
    
    pie_1=
      ggplot(data = ggplot_drug_perc[order(ggplot_drug_perc$perc,decreasing=T),],
             aes(x = "plot",
                 y = perc)
      ) +
      geom_bar(
        aes(fill = reorder(drug_class,perc)
        ),
        stat="identity",
        position = "stack"
      ) +
      labs(fill = "Drug Class") +
      facet_grid(~Sample) +
      scale_fill_got(discrete = TRUE, option = "Margaery",alpha = 0.9) +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()) +
      geom_text(aes(label = perc
      ),
      position = position_stack(vjust = .5))
    
    sup_pie_1 <- pie_1 + coord_polar("y", start=0)
    
    page1 <- plot_grid(p1,sup_pie_2,align = "hv")
    
    page2 <- plot_grid(p2,sup_pie_1,align = "hv")
    
  }
  
  # creat report dirs
  
  if (identical(path_out, "") == TRUE) 
  {out_dir <- paste(path_in,"Report", sep = "/")} else {out_dir <- paste(path_out,"Report", sep = "/")}
  
  if (identical(Exp_name, "") == TRUE)
  {sep_exp = ""} else {sep_exp = "_"}
  
  out_dir1 <- paste(out_dir,Exp_name, sep = sep_exp)
  dir.create(out_dir1, showWarnings = FALSE)
  mainDir <- out_dir1
  subDir <- head(paste(temp_dat$Sample),1)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  final_path <- file.path(mainDir, subDir)
  final_path_tab <- file.path(mainDir, subDir,subDir)
  final_path_tab_tsv <- paste(paste(final_path_tab,"Full_Report",sep = "_"),"csv", sep = ".", collapse = NULL)
  
  # evaluation file
  
  dir.create(file.path(out_dir1, "Evaluation"), showWarnings = FALSE)
  eva_file <- file.path(out_dir1, "Evaluation")
  final_path_evat <- paste(paste(eva_file,subDir,sep = "/"),"csv", sep = ".", collapse = NULL)
  
 merge1 <-  merge(ggplot_pep_count,ggplot_names)
   merge2 <-  merge(merge1,ggplot_PQ_amr)
     merge3 <- left_join(x = merge2, y = aro_index_DIA[ , c("Names","drug_class")], by = "Names")

 merge3$ID <- seq.int(nrow(merge3))
 
 merge3_sep <- separate_rows(merge3,drug_class,sep=";",convert = TRUE)

 
 # Table report
 
 if(ggplot_drug_class$drug_class == "None" ) {} else {
 
 merge4_sep <- merge(merge3_sep,ggplot_drug_perc,by = "drug_class")
 
 merge4_sep$Sample.y <- NULL
 
 names(merge4_sep)[names(merge4_sep) == "Sample.x"] <- "Sample"
 names(merge4_sep)[names(merge4_sep) == "perc.x"] <- "Name_perc"
 names(merge4_sep)[names(merge4_sep) == "perc.y"] <- "Drug_perc"

 
 write.table(merge4_sep,file = final_path_evat,row.names=FALSE, na=" ", sep=";")
 }
 
 if(exists('quant_ARO_iso')) {write.table(quant_ARO_iso,file = final_path_tab_tsv,row.names=FALSE, na=" ", sep=";")} else 
 {}
   
  # saves PDFs  
 
 pdf_Main <- head(paste(temp_dat$Sample,"main_page.pdf",sep = "_"),1)
 ggsave(pdf_Main, plot = main_page, device = NULL, path = final_path, dpi = 300,width = 10, height = 10, units = "in")
 
  pdf_1 <- head(paste(temp_dat$Sample,"overview_page.pdf",sep = "_"),1)
  ggsave(pdf_1, plot = page1, device = NULL, path = final_path, dpi = 300,width = 26.9, height = 17.4, units = "in")
  
  pdf_2 <- head(paste(temp_dat$Sample,"PQ_top3.pdf",sep = "_"),1)
  ggsave(pdf_2, plot = page2, device = NULL, path = final_path, dpi = 300,width = 26.9, height = 17.4, units = "in")
  
  # creates combined pdf
  pdf_path <- paste(final_path,"Report.pdf",sep = "_")
  
  MyPlots = list(main_page,page1,page2)
  
  dev.list()
  dev.cur()
  pdf(width=20,height=13,pdf_path)
  print(MyPlots)
  dev.off()
  
  # loop counter
  loop_count <- loop_count + 1
  info <- sprintf("%d%% done", round((loop_count/length(filenames))*100))
  setWinProgressBar(pb, (loop_count/length(filenames))*100, label=info)
  
  # loop end msg
  cat(underline(paste("Done with", subDir, sep=" ")) %+%
        "\n" %+%
        "\n")
  # paste end
}

# close pop up
close(pb)

# end msg
len_file <- length(filenames)
past1 <- paste("Saved a total of",len_file, set = "")
cat(bold(blue("Analyze all",paste(past1,"reports.",sep = ""), "with following options:","\n")))
cat(paste("Output in",out_dir1,sep = " "))

# clean up
while (!is.null(dev.list()))  dev.off()
rm(list=ls())