
# rank by numeric value
x <- paste0('X', 1:30)
str_sort(x, numeric = T)

# extract
str_extract('Subject 1111-2222.csv', '\\d+-\\d+')
