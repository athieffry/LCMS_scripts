##### Plotting caffeine Internal Standard areas, then output a csv with concentrations
##### November 2020

################ SCRIPT SETTINGS ##########################
# working directory                                       #
setwd('~/Dropbox/KEASLING/Lars R code/raw/')              # input files must all be in that directory
# inputs (see sections 1a and 2a)                         #
input_file_csv <- '210331_tmet496_sc.csv'                 # name of input file (must be csv format)
threshold_file_xlsx <- 'thresholds.xlsx'                  # name of threshold file (must be .xlsx format, do not change its column names, be sure there is one threshold for all type of peak!)
# outputs (see section 4b)                                #
save_csv <- TRUE                                          # whether to save a CSV result file (TRUE/FALSE)
output_file_csv <- 'Result_table_20200430.csv'            # name of output file (must end with .csv)
save_xlsx <- TRUE                                         # whether to save an XLSX result file (TRUE/FALSE)
output_file_xlsx <- 'Result_table_20200430.xlsx'          # name of output file (must end with .xlsx)
output_verification_csv <- 'verification.csv'             # name of output file (must end with .csv)
# plots (see sections 2f and 3c)                          #
save_plots <- TRUE                                        # should plots be saved? (TRUE/FALSE)
filter_plot_name <- 'F1. filtering plot.pdf'              # name for threshold plot (must end with .pdf, avoid bizarre characters)
caffeine_plot_name <- 'F2. caffeine plot.pdf'             # name for caffeine IS plot (must end with .pdf, avoid bizarre characters)
verification_plot_name <- 'F3. verification plot.pdf'     # name for verification plot (must end with .pdf, avoid bizarre characters)
###########################################################



# 0. PREPARING ENVIRONMENT ####
# -----------------------------
# required libraries
library(tidyverse)      # for easy wrangling of data and plotting (ggplot2)
library(tidylog)        # outputs statistics on what you do, like how many rows were modified etc. (you'll see)
library(RColorBrewer)   # color palettes
library(magrittr)       # for piping (%>%)
library(reshape2)       # to transforming data from wide to long format
library(readxl)         # to read excel files
library(WriteXLS)       # to output XLSX file (do: install.packages('WriteXLS') if you need to install it)
library(janitor)        # to clean the raw dataset (remove empty columns)
'%notin%' <- Negate('%in%') # custom function



# 1. PREPARE DATA ####
# --------------------
# 1a. read input file
   # find number of lines to skip
   skip_n <- read_lines(input_file_csv) %>%
             str_detect('Peaks') %>%
             which() %>%
             subtract(1)
   # read file
   mia <- read.csv(input_file_csv, skip=skip_n, sep=',', quote='"', header=T, strip.white=T,
                   stringsAsFactors=F, blank.lines.skip=T, na.strings=c('NA', '---'), dec='.') %>%
          as_tibble() %>%
          remove_empty(which=c('rows', 'cols'))

# 1b. rename columns
colnames(mia) <- c('file', 'note', 'peak_name', 'RT_min', 'area', 'amount', 'calib_amount', 'dev_pc_amount', 'sample_type')

# 1c. transform amounts to numeric (non-numeric values will be set to NA)
suppressWarnings( mia %<>% mutate('amount'=as.numeric(amount),
                                  'calib_amount'=as.numeric(calib_amount),
                                  'dev_pc_amount'=as.numeric(dev_pc_amount)) )

# 1d. change all different types of blanks to the same name
blank_names <- c('Blank ', 'Solvent Blank ', 'Blank_YPD ', 'Blank_SC ', 'Double Blank ',
               'Double_Blank_SC ', 'Double_Blank_YPD ', 'Double Blank_SC ', 'Double Blank_YPD ')

mia %<>% mutate('sample_type_tr'=ifelse(note %in% blank_names, 'Blank', sample_type))



# 2. APPLY THRESHOLD FILTERING ####
# ---------------------------------
# 2a. read threshold input file
threshold <- read_xlsx(threshold_file_xlsx) %>% as_tibble()
# 2b. add threshold to mia dataframe (by peak_name)
mia %<>% left_join(threshold, by='peak_name')
# 2c. add new column "area_filtered" with NA if area is below threshold
mia %<>% mutate('area_filtered'=ifelse(area < threshold, NA, area))
# 2d. remove the threshold column that is not needed anymore
mia %<>% select(-threshold)
# 2e. plot area and area_filtered to see difference
gg_filter_plot <- mia %>%
                  select(peak_name, sample_type_tr, area, area_filtered) %>%
                  melt(id.vars=c('peak_name', 'sample_type_tr'), value.name='area', variable.name='filtered') %>%
                  mutate('passed_filter'=ifelse(filtered=='area_filtered', 'yes', 'no')) %>%
                  ggplot(aes(x=area+1, col=sample_type_tr, lty=passed_filter)) +
                         geom_density() +
                         facet_wrap(~peak_name, scale='free_y') +
                         scale_x_log10() + theme(legend.position='bottom') +
                         scale_linetype_manual(values=c(2, 1), name='Passed filter') +
                         labs(title='Area density plots: before & after filtering', x='Area') +
                         scale_color_brewer(palette='Dark2', name='Sample type')
# 2f. save plot if asked to, otherwise just plot in RStudio
if(save_plots) { ggsave(plot=gg_filter_plot, filename=filter_plot_name, width=297, height=210, units='mm')
               } else gg_filter_plot



# 3. CAFFEINE ANALYSIS ####
# -------------------------
# 3a. extract caffeine-only data
caffeine <- mia %>% filter(peak_name=='Caffeine')
# 3b. plot both area and area_filtered to see if linear
gg_caffeine_plot <- caffeine %>%
                    select(file, area, area_filtered, sample_type_tr) %>%
                    melt(id.vars=c('file', 'sample_type_tr'), variable.name='filtered') %>%
                    mutate('passed_filter'=ifelse(filtered=='area_filtered', 'Filtered', 'Raw')) %>%
                    ggplot(aes(x=file, y=value, col=sample_type_tr)) +
                           geom_point() + facet_wrap(~passed_filter, ncol=1) + theme_bw() +
                           theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
                           labs(x='Files', y='Area', title='Caffeine IS: areas before & after filtering') +
                           scale_color_brewer(palette='Dark2', name='Sample type')
# 3c. save plot if asked to, otherwise just plot in RStudio
if(save_plots) { ggsave(plot=gg_caffeine_plot, filename=caffeine_plot_name, width=297, height=210, units='mm')
               } else gg_caffeine_plot




# 4. GET CSV OUT ####
# -------------------
# 4a. remove the Internal Standard (caffeine) and keep only Analysis-type results
results <- mia %>% filter(peak_name != 'Caffeine', sample_type=='Analysis')
# 4b. write to output file
if(save_csv) { write.csv(results, file=output_file_csv, append=F, quote=F, sep=',', dec='.', row.names=F, col.names=T) }
if(save_xlsx) { WriteXLS::WriteXLS(results, ExcelFileName=output_file_xlsx, SheetNames='Results', row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T) }



# 5c. Remove caffeine and rename amount column
results %>% filter(peak_name %notin% "Caffeine") %>% rename("Amt_ug_L"='Amt_filtered')

# 5d. make a wide table
mia_result_table_wide <- mia_result_table %>% pivot_wider(names_from = peak_name, values_from = c(Amt_ug_L,Amt_uM))
results <- mia_result_table_wide







# 3. VERIFICATION CONCENTRATION PLOT ####
# ---------------------------------------
# 3a. Remove Caffeine data (IS) and retain only the "Verification" samples
verification <- mia %>% filter(peak_name!='Caffeine', sample_type_tr=='Verification')
# 3b. select relevant columns
verification %<>% select(file, note, amount, peak_name)
# 3c. plotting of the verification samples to spot trends
gg_verfication_plot <- ggplot(verification, aes(x=note, y=amount, color=peak_name)) +
                              geom_point() +
                              geom_line(aes(group=1)) +
                              facet_wrap (~peak_name, scale='free') +
                              labs(y='Analyte amount_ug/L', x='Injection') +
                              theme(axis.text.x=element_text(angle=90, size=4))
gg_verfication_plot
# 3d. save plot if asked to, otherwise just plot in RStudio
if(save_plots) { ggsave(plot=gg_verfication_plot, filename=verification_plot_name, width=297, height=210, units='mm')
               } else gg_verfication_plot

# 3e. Make a summary of the verification samples by calculating the mean, sd and RSD%
verification_summary <- verification %>%
                        group_by(peak_name) %>%
                        summarize('amount_mean'=mean(amount, na.rm=T),
                                  'amount_sd'=sd(amount, na.rm=T)) %>%
                        mutate('RSD'=amount_sd*100/amount_mean)
verification_summary

# 3f. Save the output
if(save_csv) { write.csv(verification_summary, file=output_verification_csv, append=F, quote=F, sep=',', dec='.', row.names=F, col.names=T) }



# 4. PREPARING A RESULT TABLE
# -------------------------
# a. Select only the sample and relevant columns
mia_result_table <- mia %>%
                    filter(sample_type_tr=='Analysis') %>%
                    select(file, note, Amt_filtered, peak_name)

# b. make dataframe with molecular weights
MWs <- tibble('peak_name'=c('Caffeine', 'Catharanthine', 'Loganic acid', 'Loganin', 'Rauwolscine', 'Secologanin',
                            'Serpentine', 'Strictosidine', 'Tabersonine', 'Tetrahydroalstonine', 'Tryptamine',
                            'Tryptophan', 'Vinblastine', 'Vindoline'),
       'MW'=c(194.19, 336.43, 376.36, 390.38, 354.4, 388.37, 349.4, 530.57, 336.43, 352.4, 160.22, 204.23, 811, 456.54))
MWs
# c. add MW to mia dataframe and compute uM


mia_result_table%<>%
  mutate('Amt_uM'=ifelse(peak_name == 'Caffeine' , Amt_filtered/194.19, Amt_filtered),
         'Amt_uM'=ifelse(peak_name == 'Catharanthine' , Amt_filtered/336.43, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Loganic acid' , Amt_filtered/376.36, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Loganin' , Amt_filtered/390.38, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Rauwolscine' , Amt_filtered/354.4, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Secologanin' , Amt_filtered/388.37, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Serpentine' , Amt_filtered/349.4, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Strictosidine' , Amt_filtered/530.57, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Tabersonine' , Amt_filtered/336.43, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Tetrahydroalstonine' , Amt_filtered/352.4, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Tryptamine' , Amt_filtered/160.22, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Tryptophan' , Amt_filtered/204.23, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Vinblastine' , Amt_filtered/811, Amt_uM),
         'Amt_uM'=ifelse(peak_name == 'Vindoline' , Amt_filtered/456.54, Amt_uM))

# c. Remove caffeine and rename amount column
mia_result_table%<>%filter(peak_name %notin% "Caffeine")%>%rename("Amt_ug_L"='Amt_filtered')

# d. make a wide table
mia_result_table_wide<-mia_result_table%>% pivot_wider(names_from = peak_name, values_from = c(Amt_ug_L,Amt_uM))
results<-mia_result_table_wide

if(save_csv) { write.csv(results, file=output_file_csv, append=F, quote=F, sep=',', dec='.', row.names=F, col.names=T) }
if(save_xlsx) { WriteXLS::WriteXLS(results, ExcelFileName=output_file_xlsx, SheetNames='Results', row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T) }
