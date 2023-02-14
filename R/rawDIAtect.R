#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
while (!is.null(dev.list()))  dev.off()
rm(list=ls(all=TRUE))

rawDIAtect <- function(path_in = "",
                       pep_filter = 3,
                       iso_diff = 0,
                       path_out = "",
                       Exp_name = "")
                         {
# Set path to aro index database (*.tsv)

path_db_in <- system.file("data", "aro_index_DIA.tsv", package = "rawDIAtect")
options(readr.num_columns = 0)
options(warn=-1)

  if (path_out == "") {path_out = ""} else {}
  if (Exp_name == "") {Exp_name = ""} else {}
  if (pep_filter == "") {pep_filter = 3} else {}
  if (iso_diff == "") {iso_diff = 0} else {}

#rawSplit################################################################################################################
cat("\n")

cat("\n")
cat(crayon::bold(crayon::underline("Loading module rawSPLIT")))
cat("\n")
cat("\n")

# gets filenames_tsv from path_in
filenames_tsv <- list.files(path=path_in,
                            pattern=".tsv")

l1 = length(filenames_tsv)
l2 = paste("Found",l1,sep =" ")
l3 = paste(l2,"DIA-NN Quants in")
cat(crayon::blue(paste(l3, path_in,sep=" ")))
cat("\n")
cat("\n")

if (identical(Exp_name, "") == TRUE)
{sep_exp = ""} else {sep_exp = "_"}

for(i in filenames_tsv)
{
  filepath <- file.path(path_in,paste(i,sep=""))
  t1 = paste("Loading DIA-NN", i,sep = " ")
  t2 = paste(t1, "from",sep = " ")
  print(paste(t2,path_in,sep = " "))
  temp_dat <- suppressMessages(readr::read_delim(filepath, "\t", escape_double = FALSE, trim_ws = TRUE))
  cat(crayon::bgGreen$black$bold(paste("Starts splitting", i, sep= " ")))
        cat("\n")

  #read quant
  quant <- temp_dat


  #Spliting
  Sample_raw <- basename(as.character(quant[["File.Name"]]))
  Sample <- tm::removeWords(Sample_raw, ".raw")
  Samples <- as.data.frame(Sample)
  Samples[["rownumber"]] = 1:nrow(Samples)
  quant[["rownumber"]] = 1:nrow(quant)
  quant_split <- merge(quant,Samples, by="rownumber")
  quant_split <- quant_split %>% dplyr::select(-rownumber)

  # number of Files

  file_count <- length(unique(quant_split[["Sample"]]))

  file_print <- paste(file_count,"splitted Quant files saved to", sep = " ")

  #Create Dir

  out_dir <- paste(path_in,"split/", sep = "/")
  dir.create(path_out, showWarnings = FALSE)
  dir.create(out_dir, showWarnings = FALSE)

  # Split by variable
  spt2 <- split(quant_split, quant_split[["Sample"]])


  pb <- progress::progress_bar$new(format = " progress [:bar] :percent eta: :eta", # add
                         total = length(names(spt2)), clear = FALSE, width= 60)

  # check if output file are present
  list_new <- names(spt2)
  list_new_csv=c(paste0(list_new, ".csv"))
  list_old_csv <- list.files(path=out_dir,
                             pattern=".csv")
  matches <- list_new_csv %in% list_old_csv
  df_matches <- data.frame(list_new,matches)
  names(df_matches)[1] <- "names"
  df_matches_filt <- as.data.frame(filter(df_matches, matches == "FALSE"))
  list_needed_csv <- as.character(df_matches_filt$names)

  # Save
  lapply(list_needed_csv, function(x){
    pb$tick()
    write_csv(spt2[[x]], path = paste(out_dir,x,paste( Exp_name, "quant.csv",sep ="_"), sep = ""))
  })

  # loop end msg
  cat(crayon::underline(paste("Done splitting with", i, sep=" ")))
                          cat("\n")
                        paste(file_print,out_dir,sep = " ")
                        cat("\n")
  cat("\n")
  cat("\n")

}

cat(paste("Saved all splitted Quant files in", out_dir,sep = " "))
      cat("\n" )

#rawCRIT############################################################################

cat("\n")
cat(crayon::bold(crayon::underline("Loading module CRIT.r")))
cat("\n")
cat("\n")

# gets filenames from path_in
path_split <- paste(path_in,"split",sep = "/")
filenames <- list.files(path=path_split,
                        pattern=".csv")

cat(crayon::green(paste("Found splitted Quants in", path_in,sep = " ")))
    cat("\n")
    cat(crayon::green(paste("Found Index in", path_db_in,sep = " ")))
    cat("\n")

# set to clear theme
theme_set(theme_cowplot())

# reads ARO-Index
aro_index_DIA <- suppressMessages(readr::read_delim(path_db_in, "\t", escape_double = FALSE, trim_ws = TRUE))

cat(crayon::underline$green("ARO Index was successfully loaded!"))
      cat("\n")
      cat("\n")
      cat("\n")


# starts popup progress
pb <- tcltk::tkProgressBar(title="Total Progress of rawDetect", label="0% done", min=0, max=100, initial=0)
loop_count <- 0


cat((crayon::bold$bgCyan("Starting analyzing with following options:","\n")))
cat(crayon::bold(paste("Unique Peptide filter = ", pep_filter,"\n")))
cat(crayon::bold(paste("Protein-Isoform difference threshold = ", iso_diff)),"\n","\n","\n")

# main loops start
for(i in filenames){

  # get file path from file.names path.in
  filepath <- file.path(path_split,paste(i,sep=""))

  cat(crayon::bold(paste("Loading",i,sep = " "))%+%
        "\n")

  # input from loop (filepath)
  temp_dat <-  suppressMessages(readr::read_delim(filepath, ",", escape_double = FALSE, trim_ws = TRUE))

  # paste start

  cat(crayon::bgGreen$black$bold(head(paste("Starts analyzing",head(paste(temp_dat[["Sample"]]),1),sep = " "),1)))
        cat("\n")
          cat("\n")

  # greps only ARO tags from temp
  quant_ARO <- temp_dat[grep("ARO:", temp_dat[["Genes"]]), ]

  # only dplyr::selects (Genes, Stripped.Sequence, Precursor.Quantity)
  quant_ARO_small <- dplyr::select(quant_ARO,c("Sample","Genes","Stripped.Sequence","Precursor.Quantity"))

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

  if(max(quant_ARO_count$peptides) <= pep_filter |  nrow(quant_ARO_count) == 0) {
    mock_data <- data.frame(matrix(ncol = 9, nrow = 0))
    col_name_mock <- c("Sample","Names","Precursor.Quantity.Top3","amr_gene_family","Subfamily","peptides","resistance_mechanism","drug_class","Note")
    colnames(mock_data) <- col_name_mock
    mock_data[nrow(mock_data) + 1,] = c("None",head(paste(temp_dat[["Sample"]]),1),"None",0,0,"None","None","None","None")
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
    quant_ARO_index_top3_sort <- dplyr::select(quant_ARO_index_top3,c("Sample","Stripped.Sequence","Precursor.Quantity.y","amr_gene_family","peptides"))
    names(quant_ARO_index_top3_sort)[3] <- "Precursor.Quantity.Top3"

    if (pep_filter > 1)
    {
      quant_ARO_index_top3_filt <- quant_ARO_index_top3_sort %>% filter( quant_ARO_index_top3_sort[["peptides"]] > pep_filter)
    } else {
      quant_ARO_index_top3_filt <- quant_ARO_index_top3_sort
    }

    # gets iso_names
    quant_ARO_iso_index <- merge(x = quant_ARO_sep, y = aro_index_DIA[ , c("Genes","Names","amr_gene_family","Subfamily","Note")], by = "Genes", all.x=TRUE)
    quant_ARO_iso_count <- quant_ARO_iso_index %>% group_by(amr_gene_family) %>% add_count(Names,name = "count_name")

    max_isos <- quant_ARO_iso_count %>%
                  group_by(amr_gene_family) %>%
      summarise(max_count_name = max(count_name))

    #names(max_isos)[1] <- "amr_gene_family"
    #names(max_isos)[2] <- "max_count_name"
    quant_ARO_iso_count_max <- merge(quant_ARO_iso_count,max_isos)
    quant_ARO_iso_count_diff <- quant_ARO_iso_count_max %>% mutate(iso_count_diff = max_count_name-count_name)

    iso_name_best <- quant_ARO_iso_count_diff %>% filter(as.numeric(iso_count_diff) <= iso_diff)


    #iso_name_best_small <- dplyr::select(iso_name_best,c("amr_gene_family","Genes","Names","Subfamily","Note"))
    iso_name_best_small_uni <-  unique(iso_name_best)

    # output with corrected isoform
    quant_ARO_name <- merge(quant_ARO_index_top3_filt,iso_name_best_small_uni)
    quant_ARO_iso <- left_join(x = quant_ARO_name, y = aro_index_DIA[ , c("Genes","resistance_mechanism","drug_class")], by = "Genes")
    ggplot_filter_small_drugs <- as.data.frame(quant_ARO_iso)

  }

  # if no pattern found print mock found for visulication
  if(nrow(ggplot_filter_small_drugs) == 0) {
    ggplot_filter_small_drugs[nrow(ggplot_filter_small_drugs)+1,] <- NA
    ggplot_input <- ggplot_filter_small_drugs
  } else {
    ggplot_input <- ggplot_filter_small_drugs
  }

  #data
  #ggplot_input = visulize
  #quant_ARO_iso = table_report

  # visulize data

  # predefine colors

  group.fills <-
    c(
      "16S rRNA methyltransferase (A1408)" = "#792427FF",
      "16S rRNA methyltransferase (G1405)" = "#782427FF",
      "A.baumanii mutant Lpx [Colisitn]" = "#772528FF",
      "AAC(1)" = "#772529FF",
      "AAC(2')" = "#762629FF",
      "AAC(3)" = "#76272AFF",
      "AAC(3);AAC(6')" = "#75272BFF",
      "AAC(6')" = "#75282BFF",
      "AAC(6');ANT(3'')" = "#74282CFF",
      "ABC-F ATP-binding cassette ribosomal protection protein" = "#74292DFF",
      "ACC beta-lactamase" = "#732A2DFF",
      "ACI beta-lactamase" = "#732A2EFF",
      "ACT beta-lactamase" = "#722B2FFF",
      "ADC beta-lactamase w/ carbapenamase" = "#722B2FFF",
      "ADC beta-lactamase w/o carbapenamase" = "#712C30FF",
      "AER beta-lactamase" = "#712D31FF",
      "AIM beta-lactamase" = "#702D31FF",
      "aminocoumarin resistant gyrB" = "#702E32FF",
      "aminocoumarin resistant parE" = "#6F2F33FF",
      "aminocoumarin resistant parY" = "#6F2F33FF",
      "aminosalicylate resistant dihydrofolate reductase" = "#6E3034FF",
      "aminosalicylate resistant dihydrofolate synthase" = "#6E3035FF",
      "aminosalicylate resistant thymidylate synthase" = "#6D3135FF",
      "amp acetyltransferase" = "#6D3236FF",
      "ampC-type beta-lactamase" = "#6C3237FF",
      "ANT(2'')" = "#6C3338FF",
      "ANT(3'')" = "#6B3338FF",
      "ANT(4')" = "#6B3439FF",
      "ANT(6)" = "#6A353AFF",
      "ANT(9)" = "#6A353AFF",
      "APH(2'')" = "#69363BFF",
      "APH(2'');AAC(6')" = "#69363CFF",
      "APH(3')" = "#68373CFF",
      "APH(3'')" = "#68383DFF",
      "APH(4)" = "#67383EFF",
      "APH(6)" = "#67393EFF",
      "APH(7'')" = "#663A3FFF",
      "APH(9)" = "#653A40FF",
      "AQU beta-lactamase" = "#653B40FF",
      "ARL beta-lactamase" = "#643B41FF",
      "AST beta-lactamase" = "#643C42FF",
      "ATP-binding cassette (ABC) efflux pump" = "#633D42FF",
      "ATP-binding cassette (ABC) efflux pump;major facilitator superfamily (MFS) efflux pump" = "#633D43FF",
      "ATP-binding cassette (ABC) efflux pump;major facilitator superfamily (MFS) efflux pump;resistance-nodulation-cell division (RND) efflux pump" = "#623E44FF",
      "ATP-binding cassette (ABC) efflux pump;major facilitator superfamily (MFS) efflux pump;resistance-nodulation-cell division (RND) efflux pump;General Bacterial Porin w/ reduced permeability to beta-lactams" = "#623E44FF",
      "ATP-binding cassette (ABC) efflux pump;pmr phosphoethanolamine transferase" = "#613F45FF",
      "Bah amidohydrolase" = "#614046FF",
      "BAT beta-lactamase" = "#604046FF",
      "Bc beta-lactamase" = "#604147FF",
      "BCL beta-lactamase" = "#5F4148FF",
      "BEL beta-lactamase" = "#5F4249FF",
      "BES beta-lactamase" = "#5E4349FF",
      "BIC beta-lactamase" = "#5E434AFF",
      "BIL beta-lactamase" = "#5D444BFF",
      "BJP beta-lactamase" = "#5D454BFF",
      "BKC beta-lactamase" = "#5C454CFF",
      "BlaA beta-lactamase" = "#5C464DFF",
      "BlaB beta-lactamase" = "#5B464DFF",
      "blaF family beta-lactamase" = "#5B474EFF",
      "blaZ beta-lactamase" = "#5A484FFF",
      "Bleomycin resistant protein" = "#5A484FFF",
      "BPU beta-lactamase" = "#594950FF",
      "BRO beta-lactamase" = "#594951FF",
      "BUT beta-lactamase" = "#584A51FF",
      "CAM beta-lactamase" = "#584B52FF",
      "CARB beta-lactamase" = "#574B53FF",
      "CarO porin" = "#574C53FF",
      "CAU beta-lactamase" = "#564C54FF",
      "CblA beta-lactamase" = "#564D55FF",
      "CBP beta-lactamase" = "#554E55FF",
      "CcrA beta-lactamase" = "#554E56FF",
      "cepA beta-lactamase" = "#544F57FF",
      "CepS beta-lactamase" = "#545058FF",
      "Cfr 23S ribosomal RNA methyltransferase" = "#535058FF",
      "CfxA beta-lactamase" = "#525159FF",
      "CGA beta-lactamase" = "#52525AFF",
      "CGB beta-lactamase" = "#51525AFF",
      "chloramphenicol acetyltransferase (CAT)" = "#51535BFF",
      "chloramphenicol phosphotransferase" = "#50545CFF",
      "CIA beta-lactamase" = "#50545DFF",
      "ciprofloxacin phosphotransferase" = "#4F555DFF",
      "CKO beta-lactamase" = "#4E565EFF",
      "class A Bacillus anthracis Bla beta-lactamase" = "#4E565FFF",
      "class A LRA beta-lactamase" = "#4D5760FF",
      "class C LRA beta-lactamase" = "#4D5861FF",
      "class C LRA beta-lactamase;class D LRA beta-lactamase" = "#4C5861FF",
      "CME beta-lactamase" = "#4C5962FF",
      "CMH beta-lactamase" = "#4B5A63FF",
      "CMY beta-lactamase" = "#4A5A63FF",
      "cpa acetyltransferase" = "#4A5B64FF",
      "CphA beta-lactamase" = "#495C65FF",
      "CPS beta-lactamase" = "#495C66FF",
      "CTX-M beta-lactamase" = "#485D66FF",
      "cya adenylate cyclase" = "#485E67FF",
      "dapsone resistant dihydropteroate synthase folP" = "#475E68FF",
      "daptomycin resistant agrA" = "#465F69FF",
      "daptomycin resistant beta-subunit of RNA polymerase (rpoB)" = "#46606AFF",
      "daptomycin resistant CdsA" = "#45616AFF",
      "daptomycin resistant cls" = "#45616BFF",
      "daptomycin resistant drmA" = "#44626CFF",
      "daptomycin resistant gdpD" = "#44636CFF",
      "daptomycin resistant gshF" = "#43636DFF",
      "daptomycin resistant liaF" = "#42646EFF",
      "daptomycin resistant liaR" = "#42656FFF",
      "daptomycin resistant liaS" = "#416570FF",
      "daptomycin resistant mprF" = "#416670FF",
      "daptomycin resistant pgsA" = "#406771FF",
      "daptomycin resistant walK" = "#406772FF",
      "daptomycin resistant YvlB" = "#3F6873FF",
      "daptomycin resistant YybT" = "#3E6973FF",
      "daptomycin-resistant beta prime subunit of RNA polymerase (rpoC)" = "#3E6974FF",
      "defensin resistant mprF" = "#3D6A75FF",
      "DES beta-lactamase" = "#3D6B75FF",
      "DHA beta-lactamase" = "#3C6B76FF",
      "DIM beta-lactamase" = "#3C6C77FF",
      "DnaA chromosomal replication initiation protein" = "#3B6D78FF",
      "EBR beta-lactamase" = "#3A6D79FF",
      "Edeine acetyltransferase" = "#3A6E79FF",
      "elfamycin resistant EF-Tu" = "#396F7AFF",
      "Erm 23S ribosomal RNA methyltransferase" = "#396F7BFF",
      "ERP beta-lactamase" = "#38707CFF",
      "ESP beta-lactamase" = "#38717CFF",
      "ethambutol resistant arabinosyltransferase" = "#37727DFF",
      "Ethambutol resistant iniA" = "#36727EFF",
      "Ethambutol resistant iniB" = "#36737FFF",
      "Ethambutol resistant iniC" = "#35747FFF",
      "ethionamide resistant ethA" = "#357480FF",
      "EXO beta-lactamase" = "#347581FF",
      "FAR beta-lactamase" = "#347682FF",
      "FEZ beta-lactamase" = "#337682FF",
      "FIM beta-lactamase" = "#327783FF",
      "fluoroquinolone resistant gyrA" = "#327884FF",
      "fluoroquinolone resistant gyrB" = "#317885FF",
      "fluoroquinolone resistant parC" = "#317985FF",
      "fluoroquinolone resistant parE" = "#307A86FF",
      "Fom phosphotransferase family" = "#307A87FF",
      "FONA beta-lactamase" = "#2F7B88FF",
      "fosC phosphotransferase family" = "#2E7C88FF",
      "fosfomycin thiol transferase" = "#2E7C89FF",
      "FOX beta-lactamase" = "#2D7D8AFF",
      "FPH beta-lactamase" = "#2D7E8BFF",
      "FRI beta-lactamase" = "#2C7E8BFF",
      "FTU beta-lactamase" = "#2C7F8CFF",
      "fusidic acid inactivation enzyme" = "#2B808DFF",
      "General Bacterial Porin w/ reduced permeability to beta-lactams" = "#2B818EFF",
      "General Bacterial Porin w/ reduced permeability to peptides" = "#2C818EFF",
      "GES beta-lactamase" = "#2D818EFF",
      "GIM beta-lactamase" = "#2E828EFF",
      "gimA family macrolide glycosyltransferase" = "#2F828EFF",
      "GlpT" = "#30838EFF",
      "GOB beta-lactamase" = "#32838EFF",
      "HERA beta-lactamase" = "#33848EFF",
      "HMB beta-lactamase" = "#34848FFF",
      "IMI beta-lactamase" = "#35848FFF",
      "IMP beta-lactamase" = "#36858FFF",
      "IND beta-lactamase" = "#37858FFF",
      "intrinsic colistin resistant phosphoethanolamine transferase" = "#39868FFF",
      "isoniazid resistant katG" = "#3A868FFF",
      "JOHN beta-lactamase" = "#3B878FFF",
      "kdpDE" = "#3C8790FF",
      "KHM beta-latamase" = "#3D8790FF",
      "KPC beta-lactamase" = "#3F8890FF",
      "L1 family beta-lactamase" = "#408890FF",
      "LAP beta-lactamase" = "#418990FF",
      "LCR beta-lactamase" = "#428990FF",
      "LEN beta-lactamase" = "#438A90FF",
      "lincosamide nucleotidyltransferase (LNU)" = "#448A91FF",
      "lipid A acyltransferase" = "#468A91FF",
      "lipid A phosphatase" = "#478B91FF",
      "Llm 23S ribosomal RNA methyltransferase" = "#488B91FF",
      "lysocin resistant menA" = "#498C91FF",
      "macrolide esterase" = "#4A8C91FF",
      "macrolide phosphotransferase (MPH)" = "#4C8D91FF",
      "major facilitator superfamily (MFS) efflux pump" = "#4D8D92FF",
      "major facilitator superfamily (MFS) efflux pump;resistance-nodulation-cell division (RND) efflux pump" = "#4E8D92FF",
      "MCR phosphoethanolamine transferase" = "#4F8E92FF",
      "methicillin resistant PBP2" = "#508E92FF",
      "mgt macrolide glycotransferase" = "#518F92FF",
      "MipA-interacting Protein" = "#538F92FF",
      "MIR beta-lactamase" = "#549092FF",
      "MOX beta-lactamase" = "#559092FF",
      "MSI beta-lactamase" = "#569093FF",
      "MSI-OXA family beta-lactamase" = "#579193FF",
      "multidrug and toxic compound extrusion (MATE) transporter" = "#599193FF",
      "murA transferase" = "#5A9293FF",
      "murG transferase" = "#5B9293FF",
      "MUS beta-lactamase" = "#5C9393FF",
      "NDM beta-lactamase" = "#5D9393FF",
      "NmcA beta-lactamase" = "#5E9394FF",
      "non-erm 23S ribosomal RNA methyltransferase (A1067)" = "#609494FF",
      "non-erm 23S ribosomal RNA methyltransferase (G748)" = "#619494FF",
      "NPS beta-lactamase" = "#629594FF",
      "OCH beta-lactamase" = "#639594FF",
      "OKP beta-lactamase" = "#649694FF",
      "ole glycosyltransferase" = "#669694FF",
      "Outer Membrane Porin (Opr)" = "#679695FF",
      "OXA beta-lactamase" = "#689795FF",
      "OXY beta-lactamase" = "#699795FF",
      "PDC beta-lactamase" = "#6A9895FF",
      "Penicillin-binding protein mutations conferring resistance to beta-lactams" = "#6B9895FF",
      "PER beta-lactamase" = "#6D9995FF",
      "phenotypic variant regulator" = "#6E9995FF",
      "pmr phosphoethanolamine transferase" = "#6F9996FF",
      "PNGM beta-lactamase" = "#709A96FF",
      "PtsI phosphotransferase" = "#719A96FF",
      "pyrazinamide resistant pncA" = "#739B96FF",
      "pyrazinamide resistant rpsA" = "#749B96FF",
      "quinolone resistance protein (qnr)" = "#759C96FF",
      "R39 beta-lactamase" = "#769C96FF",
      "RbpA bacterial RNA polymerase-binding protein" = "#779C97FF",
      "RCP beta-lactamase" = "#789D97FF",
      "resistance fabG" = "#7A9D97FF",
      "resistance-nodulation-cell division (RND) efflux pump" = "#7B9E97FF",
      "resistance-nodulation-cell division (RND) efflux pump;General Bacterial Porin w/ reduced permeability to beta-lactams" = "#7C9E97FF",
      "resistant fabI" = "#7D9F97FF",
      "resistant fusA" = "#7E9F97FF",
      "resistant fusE" = "#80A098FF",
      "resistant gidB" = "#81A098FF",
      "resistant inhA" = "#82A098FF",
      "resistant isoleucyl-tRNA synthetase (ileS)" = "#83A198FF",
      "resistant kasA" = "#84A198FF",
      "resistant ndh" = "#85A298FF",
      "resistant nfsA" = "#86A298FF",
      "resistant nfsB" = "#87A298FF",
      "resistant rpsL" = "#88A399FF",
      "resistant tlyA" = "#8AA399FF",
      "rifampin ADP-ribosyltransferase (Arr)" = "#8BA499FF",
      "rifampin glycosyltransferase" = "#8CA499FF",
      "rifampin monooxygenase" = "#8DA499FF",
      "rifampin phosphotransferase" = "#8EA599FF",
      "rifamycin-resistant arabinosyltransferase" = "#8FA599FF",
      "rifamycin-resistant beta-subunit of RNA polymerase (rpoB)" = "#90A69AFF",
      "Rm3 family beta-lactamase" = "#91A69AFF",
      "ROB beta-lactamase" = "#93A69AFF",
      "RSA beta-lactamase" = "#94A79AFF",
      "Sed beta-lactamase" = "#95A79AFF",
      "SFH beta-lactamase" = "#96A89AFF",
      "SHV beta-lactamase" = "#97A89AFF",
      "SHW beta-lactamase" = "#98A89BFF",
      "SIM beta-lactamase" = "#99A99BFF",
      "small multidrug resistance (SMR) efflux pump" = "#9AA99BFF",
      "SMB beta-lactamase" = "#9CAA9BFF",
      "SME beta-lactamase" = "#9DAA9BFF",
      "SPG beta-lacatamase" = "#9EAA9BFF",
      "SPM beta-lactamase" = "#9FAB9BFF",
      "SRT beta-lactamase" = "#A0AB9CFF",
      "streptogramin vat acetyltransferase" = "#A1AC9CFF",
      "streptogramin vgb lyase" = "#A2AC9CFF",
      "streptothricin acetyltransferase (SAT)" = "#A3AC9CFF",
      "subclass B1 Bacillus anthracis Bla beta-lactamase" = "#A5AD9CFF",
      "subclass B1 PEDO beta-lactamase" = "#A6AD9CFF",
      "subclass B1 Vibrio cholerae varG beta-lactamase" = "#A7AE9CFF",
      "subclass B3 LRA beta-lactamase" = "#A8AE9CFF",
      "subclass B3 PEDO beta-lactamase" = "#A9AE9DFF",
      "Sugar Porin (SP)" = "#AAAF9DFF",
      "sulfonamide resistant dihydropteroate synthase folP" = "#ABAF9DFF",
      "sulfonamide resistant sul" = "#ACB09DFF",
      "TEM beta-lactamase" = "#AEB09DFF",
      "tetracycline inactivation enzyme" = "#AFB09DFF",
      "tetracycline-resistant ribosomal protection protein" = "#B0B19DFF",
      "THIN-B beta-lactamase" = "#B1B19EFF",
      "TLA beta-lactamase" = "#B2B29EFF",
      "TMB beta-lactamase" = "#B3B29EFF",
      "triclosan resistant gyrA" = "#B4B29EFF",
      "trimethoprim resistant dihydrofolate reductase dfr" = "#B5B39EFF",
      "TRU beta-lactamase" = "#B7B39EFF",
      "tunicamycin resistance protein" = "#B8B49EFF",
      "TUS beta-lactamase" = "#B9B49FFF",
      "UhpT" = "#BAB49FFF",
      "undecaprenyl pyrophosphate related proteins" = "#BBB59FFF",
      "van ligase" = "#BCB59FFF",
      "vancomycin-resistant beta prime subunit of RNA polymerase (rpoC)" = "#BDB69FFF",
      "vanH" = "#BEB69FFF",
      "vanJ membrane protein" = "#C0B69FFF",
      "vanK" = "#C1B7A0FF",
      "vanR" = "#C2B7A0FF",
      "vanS" = "#C3B8A0FF",
      "vanT" = "#C4B8A0FF",
      "vanU" = "#C5B8A0FF",
      "vanV" = "#C6B9A0FF",
      "vanW" = "#C8B9A0FF",
      "vanX" = "#C9BAA1FF",
      "vanXY" = "#CABAA1FF",
      "vanY" = "#CBBAA1FF",
      "vanZ" = "#CCBBA1FF",
      "VCC beta-lactamase" = "#CDBBA1FF",
      "VEB beta-lactamase" = "#CEBCA1FF",
      "VIM beta-lactamase" = "#CFBCA1FF",
      "viomycin phosphotransferase" = "#D1BDA2FF"
    )

  group.colors <- c(broadspectrum = "#006400",
                    Carbapenemase = "#DC143C",
                    ESBL ="#ffcccb",
                    Narrow = "#90ee90",
                    narrowspectrum = "#32CD32")

  # unique peptide count

  ggplot_pep_count <- unique(dplyr::select(ggplot_input,c("Sample","amr_gene_family","peptides")))
  ggplot_pep_names <- unique(dplyr::select(ggplot_input,c("Sample","amr_gene_family","Names","Subfamily","Note")))
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
    facet_grid(~head(temp_dat[["Sample"]],1)) +
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

  ggplot_PQ_amr <- unique(dplyr::select(ggplot_input,c("Sample","amr_gene_family","Precursor.Quantity.Top3")))

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
  amr_gene_family <- c("AAC(6')","AAC(3')","ANT(6)",
                       "KPC beta-lactamase","OXA beta-lactamase",
                       "SHV beta-lactamase","TEM beta-lactamase",
                       "CMY beta-lactamase","CTX-M beta-lactamase",
                       "NDM beta-lactamase","VIM beta-lactamase",
                       "MCR phosphoethanolamine transferase",
                       "van ligase","vanR","vanS")
  
  group_amr <- c("Aminoglycoside","Aminoglycoside","Aminoglycoside",
                 "Beta-lactam","Beta-lactam",
                 "Beta-lactam","Beta-lactam",
                 "Beta-lactam","Beta-lactam",
                 "Beta-lactam","Beta-lactam",
                 "Colistin",
                 "Vancomycin","Vancomycin","Vancomycin")
  Found_AMR <- c("")
  main_page_dat <- data.frame(amr_gene_family,group_amr,Found_AMR)
  
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
    ggplot2::geom_tile(aes(x = amr_gene_family,
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
    facet_grid(Note~group_amr, scales = "free") +
  labs(title = head(temp_dat[["Sample"]],1))


  # drug_class and spectrum

  ggplot_drug_class <- dplyr::select(ggplot_input,c("Sample","drug_class","peptides"))
  ggplot_drug_class <- as.tbl(unique(separate_rows(ggplot_drug_class,drug_class,sep=";",convert = TRUE)))


  # if no drug classes

  if(length(ggplot_drug_class$drug_class)== 0 | is.na(ggplot_drug_class$drug_class)) {
    page1 <- p1
    page2 <- p2

  } else {
    # gets sums and percent of drug_classes

    ggplot_drug_sum <- aggregate(as.numeric(as.character(ggplot_drug_class$peptides)), by=list(drug_class=ggplot_drug_class$drug_class), FUN=sum)
    ggplot_drug_perc <- ggplot_drug_sum %>% mutate(perc = round((x/sum(x))*100,2))
    ggplot_drug_perc[["Sample"]] <- rep(head(ggplot_drug_class[["Sample"]],1),nrow(ggplot_drug_perc)) # make new column

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
      facet_grid(~head(temp_dat[["Sample"]],1)) +
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

  out_dir1 <- paste(out_dir,Exp_name, sep = sep_exp)
  dir.create(out_dir1, showWarnings = FALSE)
  mainDir <- out_dir1
  subDir <- head(paste(temp_dat[["Sample"]]),1)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  final_path <- file.path(mainDir, subDir)
  final_path_tab <- file.path(mainDir, subDir,subDir)
  final_path_tab_tsv <- paste(paste(final_path_tab,Exp_name,"Full_Report" ,sep = "_"),"csv", sep = ".", collapse = NULL)

  if(exists('quant_ARO_iso')) {write.table(quant_ARO_iso,file = final_path_tab_tsv,row.names=FALSE, na=" ", sep=";")} else
  {}

  # saves PDFs

  # creates combined pdf
  pdf_path <- paste(final_path, paste(Exp_name, "Report.pdf",sep = "_"), sep = "_")

  MyPlots = list(main_page,page1,page2)

  dev.list()
  dev.cur()
  pdf(width=20,height=13,pdf_path)
  print(MyPlots)
  dev.off()

  # loop counter
  loop_count <- loop_count + 1
  info <- sprintf("%d%% done", round((loop_count/length(filenames))*100))
  tcltk::setTkProgressBar(pb, (loop_count/length(filenames))*100, label=info)

  # loop end msg
  cat(crayon::underline(paste("Done with", subDir, sep=" ")))
        cat( "\n")
        cat( "\n")
  # paste end
        rm(list=setdiff(ls(all=TRUE), c( "path_in", "pep_filter", "iso_diff", "path_out", "Exp_name",
                                 "path_split", "filenames",  "aro_index_DIA", "pb", "loop_count", "out_dir1", "sep_exp")))
}

# close pop up
close(pb)

# end msg
len_file <- length(filenames)
past1 <- paste("Saved a total of",len_file, set = "")
cat(crayon::bold(crayon::blue(paste(past1,"reports.",sep = ""), "Run with following options:","\n","\n")))
cat(paste("Unique Peptide filter = ", pep_filter,"\n"))
cat(paste("Protein-Isoform difference threshold = ", iso_diff,"\n"))
cat(paste("Output saved",out_dir1,sep = " "),"\n","\n")
cat(crayon::bold(crayon::underline("Please cite:","\n","\n")))
cat(paste("Unbiased Antimicrobial Resistance Detection from Clinical Bacterial Isolates Using Proteomics"),"\n")
cat(paste("Christian Blumenscheit, Yvonne Pfeifer, Guido Werner, Charlyn John, Andy Schneider, Peter Lasch, and Joerg Doellinger"),"\n")
cat(paste("Analytical Chemistry 2021 93 (44), 14599-14608"),"\n")
cat(paste("DOI: 10.1021/acs.analchem.1c00594")) 


# clean up
while (!is.null(dev.list()))  dev.off()
rm(list=ls(all=TRUE))

}
