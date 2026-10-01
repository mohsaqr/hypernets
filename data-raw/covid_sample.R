# Build the bundled `covid_sample` dataset: a simple random sample of 1,000
# abstracts from the COVID-19 education-research Scopus export also used by
# the sbert package (source file: ../SBERT/covid.csv, 4,187 records,
# 2020-2024). Records without a usable abstract are excluded first, and a
# repeated abstract is kept once.
# Run from the package root: Rscript data-raw/covid_sample.R

raw <- read.csv(file.path("..", "SBERT", "covid.csv"))

usable <- subset(
  raw,
  !is.na(Abstract) & nchar(Abstract) >= 400 &
    Abstract != "[No abstract available]" & !is.na(Year),
  select = c("EID", "Title", "Abstract", "Year")
)
usable <- usable[!duplicated(usable$Abstract), ]
# EID order before sampling, so the draw does not depend on the file order
usable <- usable[order(usable$EID), ]

set.seed(20261001)
picked <- usable[sort(sample(nrow(usable), size = 1000L)), ]

covid_sample <- data.frame(
  doc = picked$EID,
  title = picked$Title,
  abstract = picked$Abstract,
  year = as.integer(picked$Year),
  row.names = NULL
)

stopifnot(
  "1,000 abstracts" = nrow(covid_sample) == 1000L,
  "no duplicate document IDs" = anyDuplicated(covid_sample$doc) == 0L,
  "no duplicate abstracts" = anyDuplicated(covid_sample$abstract) == 0L,
  "no missing abstracts" = !anyNA(covid_sample$abstract)
)

dir.create("data", showWarnings = FALSE)
save(covid_sample, file = file.path("data", "covid_sample.rda"),
     compress = "xz")
cat("rows:", nrow(covid_sample), "| usable records:", nrow(usable),
    "| size:", file.size(file.path("data", "covid_sample.rda")), "bytes\n")
