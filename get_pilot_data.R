#### Get pilot/testing data for BDC-NSRR
#### Diego Mazzotti
#### October 2020

# Access and download 3 EDFs and annotations from each approved study (abc, bestair, ccshs, chat, hchs, heartbeat, homepap, mesa, mros, numom2b, shhs, sof). Access to msc is currently under review.

# Load packages
library(nsrr)
library(dplyr)

# Authenticate TOKEN (stored in .Renviron)
nsrr_auth()

# list all files for each study
cohorts <- c("abc", "bestair", "ccshs", "chat", "hchs", "heartbeat", "homepap", "mesa", "mros", "numom2b", "shhs", "sof")

file_list <- list()

#abc
file_list$edf$abc$baseline <- nsrr_dataset_files("abc", path = "polysomnography/edfs/baseline")
file_list$edf$abc$month09 <- nsrr_dataset_files("abc", path = "polysomnography/edfs/month09")
file_list$edf$abc$month18 <- nsrr_dataset_files("abc", path = "polysomnography/edfs/month18")
file_list$anno$abc$baseline <- nsrr_dataset_files("abc", path = "polysomnography/annotations-events-nsrr/baseline")
file_list$anno$abc$month09 <- nsrr_dataset_files("abc", path = "polysomnography/annotations-events-nsrr/month09")
file_list$anno$abc$month18 <- nsrr_dataset_files("abc", path = "polysomnography/annotations-events-nsrr/month18")

#bestair
file_list$edf$bestair$baseline <- nsrr_dataset_files("bestair", path = "polysomnography/edfs/baseline")
file_list$edf$bestair$followup <- nsrr_dataset_files("bestair", path = "polysomnography/edfs/followup")
file_list$edf$bestair$nonrandomized <- nsrr_dataset_files("bestair", path = "polysomnography/edfs/nonrandomized")
file_list$anno$bestair$baseline <- nsrr_dataset_files("bestair", path = "polysomnography/annotations-events-nsrr/baseline")
file_list$anno$bestair$followup <- nsrr_dataset_files("bestair", path = "polysomnography/annotations-events-nsrr/followup")
file_list$anno$bestair$nonrandomized <- nsrr_dataset_files("bestair", path = "polysomnography/annotations-events-nsrr/nonrandomized")

#ccshs
file_list$edf$ccshs$baseline <- nsrr_dataset_files("ccshs", path = "polysomnography/edfs")
file_list$anno$ccshs$baseline <- nsrr_dataset_files("ccshs", path = "polysomnography/annotations-events-nsrr")

#chat
file_list$edf$chat$baseline <- nsrr_dataset_files("chat", path = "polysomnography/edfs/baseline")
file_list$edf$chat$followup <- nsrr_dataset_files("chat", path = "polysomnography/edfs/followup")
file_list$edf$chat$nonrandomized <- nsrr_dataset_files("chat", path = "polysomnography/edfs/nonrandomized")
file_list$anno$chat$baseline <- nsrr_dataset_files("chat", path = "polysomnography/annotations-events-nsrr/baseline")
file_list$anno$chat$followup <- nsrr_dataset_files("chat", path = "polysomnography/annotations-events-nsrr/followup")
file_list$anno$chat$nonrandomized <- nsrr_dataset_files("chat", path = "polysomnography/annotations-events-nsrr/nonrandomized")

#hchs
file_list$edf$hchs$baseline <- nsrr_dataset_files("hchs", path = "polysomnography/edfs")
file_list$anno$hchs$baseline <- nsrr_dataset_files("hchs", path = "polysomnography/annotations-events-nsrr")

#heartbeat
file_list$edf$heartbeat$baseline <- nsrr_dataset_files("heartbeat", path = "polysomnography/edfs/baseline")
file_list$edf$heartbeat$followup <- nsrr_dataset_files("heartbeat", path = "polysomnography/edfs/followup")
file_list$anno$heartbeat$baseline <- nsrr_dataset_files("heartbeat", path = "polysomnography/annotations-events-nsrr/baseline")
file_list$anno$heartbeat$followup <- nsrr_dataset_files("heartbeat", path = "polysomnography/annotations-events-nsrr/followup")

#homepap
file_list$edf$homepap$home <- nsrr_dataset_files("homepap", path = "polysomnography/edfs/home")
file_list$edf$homepap$lab_full <- nsrr_dataset_files("homepap", path = "polysomnography/edfs/lab/full")
file_list$edf$homepap$lab_split <- nsrr_dataset_files("homepap", path = "polysomnography/edfs/lab/split")
file_list$anno$homepap$home <- nsrr_dataset_files("homepap", path = "polysomnography/annotations-events-nsrr/home")
file_list$anno$homepap$lab_full <- nsrr_dataset_files("homepap", path = "polysomnography/annotations-events-nsrr/lab/full")
file_list$anno$homepap$lab_split <- nsrr_dataset_files("homepap", path = "polysomnography/annotations-events-nsrr/lab/split")

#mesa
file_list$edf$mesa$baseline <- nsrr_dataset_files("mesa", path = "polysomnography/edfs")
file_list$anno$mesa$baseline <- nsrr_dataset_files("mesa", path = "polysomnography/annotations-events-nsrr")

#mros
file_list$edf$mros$visit1 <- nsrr_dataset_files("mros", path = "polysomnography/edfs/visit1")
file_list$edf$mros$visit2 <- nsrr_dataset_files("mros", path = "polysomnography/edfs/visit2")
file_list$anno$mros$visit1 <- nsrr_dataset_files("mros", path = "polysomnography/annotations-events-nsrr/visit1")
file_list$anno$mros$visit2 <- nsrr_dataset_files("mros", path = "polysomnography/annotations-events-nsrr/visit2")

#numom2b
file_list$edf$numom2b$visit1 <- nsrr_dataset_files("numom2b", path = "polysomnography/edfs/visit1")
file_list$edf$numom2b$visit3 <- nsrr_dataset_files("numom2b", path = "polysomnography/edfs/visit3")
file_list$anno$numom2b$visit1 <- nsrr_dataset_files("numom2b", path = "polysomnography/annotations-events-profusion/visit1")
file_list$anno$numom2b$visit3 <- nsrr_dataset_files("numom2b", path = "polysomnography/annotations-events-profusion/visit3")

#shhs
file_list$edf$shhs$shhs1 <- nsrr_dataset_files("shhs", path = "polysomnography/edfs/shhs1")
file_list$edf$shhs$shhs2 <- nsrr_dataset_files("shhs", path = "polysomnography/edfs/shhs2")
file_list$anno$shhs$shhs1 <- nsrr_dataset_files("shhs", path = "polysomnography/annotations-events-profusion/shhs1")
file_list$anno$shhs$shhs2 <- nsrr_dataset_files("shhs", path = "polysomnography/annotations-events-profusion/shhs2")

#sof
file_list$edf$sof$baseline <- nsrr_dataset_files("sof", path = "polysomnography/edfs")
file_list$anno$sof$baseline <- nsrr_dataset_files("sof", path = "polysomnography/annotations-events-nsrr")


### Combine in one data frame (edfs and annotations separately)
files_by_study <- list()
for (s in cohorts) {
        
        files_by_study$edf[[s]] <- bind_rows(file_list$edf[[s]], .id = "wave")
        files_by_study$anno[[s]] <- bind_rows(file_list$anno[[s]], .id = "wave")
        
}

nsrr.edf.df <- bind_rows(files_by_study$edf, .id = "study")
nsrr.anno.df <- bind_rows(files_by_study$anno, .id = "study")


# Select a random samples of 3 EDFs and corresponding annotations for each study and wave
pilot_list <- list()
for (s in cohorts) {
        
        unique_waves <- nsrr.edf.df %>%
                filter(dataset==s) %>%
                select(wave) %>%
                unique() %>%
                pull()
        
        
        for (w in unique_waves) {
                
                done=F
                
                while (!done) {
                        
                        #set.seed(10122020)
                        edf.s <- nsrr.edf.df %>%
                                filter(dataset==s, wave==w) %>%
                                select(file_name) %>%
                                sample_n(size = 3) %>%
                                pull() %>%
                                sort()
                        
                        prefix <- sapply(strsplit(edf.s, "\\."), "[[", 1)
                        
                        # Find if xml is available
                        xml.all <- nsrr.anno.df %>%
                                filter(dataset==s) %>%#, wave==w) %>%
                                select(file_name) %>%
                                pull()
                        
                        xml.s <- grep(paste(prefix, collapse = "|"), nsrr.anno.df$file_name, value = T)
                        
                        if(length(xml.s)==3) {
                                done=T
                        }

                }
                
                
                pilot_list[[s]][[w]][["edf"]] <- edf.s
                pilot_list[[s]][[w]][["anno"]] <- xml.s
                
                
                
        }
        
        
}




edf.files.sample <- grep(".edf", unlist(pilot_list, use.names = F), value = T)
anno.files.sample <- grep(".xml", unlist(pilot_list, use.names = F), value = T)

sample.edfs <- filter(nsrr.edf.df, file_name %in% edf.files.sample)
sample.anno <- filter(nsrr.anno.df, file_name %in% anno.files.sample)



### Save Pilot Test File lists
saveRDS(nsrr.edf.df,"all_nsrr.edf.df.Rdata")
saveRDS(nsrr.anno.df,"all_nsrr.edf.df.Rdata")
saveRDS(sample.edfs,"sample.edfs.Rdata")
saveRDS(sample.anno,"sample.anno.Rdata")


#### Access data ALL DATA - not RUN
# for (f in sample.edfs$full_path) {
#         
#         
#         fname <- sample.edfs$file_name[sample.edfs$full_path==f]
#         dset <- sapply(strsplit(fname, "-"), "[[", 1)
#         
#         
#         dl_file <- nsrr_download_file(dset, path = f)
#         file.copy(from = dl_file$outfile, to = paste0("pilot_data/edf/", fname))
#         
#         file.remove(dl_file$outfile)
# }
# 
# for (f in sample.anno$full_path) {
#         
#         
#         fname <- sample.anno$file_name[sample.anno$full_path==f]
#         dset <- sapply(strsplit(fname, "-"), "[[", 1)
#         
#         
#         dl_file <- nsrr_download_file(dset, path = f)
#         file.copy(from = dl_file$outfile, to = paste0("pilot_data/annp/", fname))
#         
#         file.remove(dl_file$outfile)
# }
# 

# Select one EDF per study
sample.sample.ids <- c(1,10, 19, 23, 31, 34, 40, 49, 52, 58, 64, 70)
for (f in sample.edfs$full_path[sample.sample.ids]) {


        fname <- sample.edfs$file_name[sample.edfs$full_path==f]
        dset <- sapply(strsplit(fname, "-"), "[[", 1)
        if(dset=="shhs1") {dset <- "shhs"}


        dl_file <- nsrr_download_file(dset, path = f)
        file.copy(from = dl_file$outfile, to = paste0("pilot_data/edf/", fname))

        file.remove(dl_file$outfile)
}

for (f in sample.anno$full_path[sample.sample.ids]) {


        fname <- sample.anno$file_name[sample.anno$full_path==f]
        dset <- sapply(strsplit(fname, "-"), "[[", 1)
        if(dset=="shhs1") {dset <- "shhs"}


        dl_file <- nsrr_download_file(dset, path = f)
        file.copy(from = dl_file$outfile, to = paste0("pilot_data/anno/", fname))

        file.remove(dl_file$outfile)
}


