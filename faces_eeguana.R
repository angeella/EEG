library(httr)
GET("https://osf.io/wv5hp//?action=download",
    write_disk("./s1_faces.vhdr", overwrite = TRUE),
    progress()
)
GET("https://osf.io/c8wea//?action=download",
    write_disk("./s1_faces.vmrk", overwrite = TRUE),
    progress()
)
GET("https://osf.io/uhsde//?action=download",
    write_disk("./s1_faces.eeg", overwrite = TRUE),
    progress()
)
faces1 <- read_vhdr("faces.vhdr")

events_tbl(faces)
channels_tbl(faces)
segments_tbl(faces)

faces_seg <- faces %>%
    eeg_segment(.description %in% c("s70","s71"),
                lim = c(-0.2, 0.25)) 

