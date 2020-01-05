
fls <- list.files('data/rData/')
idx <- grep('Sample_Information', fls)

for (fl in fls[idx]) {
 file.rename(file.path('data/rData/',fl), file.path('data/rData/',gsub('Sample_Information', 'Metadata', fl)))
}
