# Build the Legal Hypergraphs datasets (Coupette, Hartung & Katz 2024) from
# the authors' Zenodo replication archive, record 8081507 (CC BY-NC 4.0).
#   Rscript data-raw/legal_hypergraphs.R /path/to/unpacked/archive
archive <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(dir.exists(file.path(archive, "data_new")))
read <- function(file, ...) {
  utils::read.csv(file.path(archive, "data_new", file), stringsAsFactors = FALSE,
                  encoding = "UTF-8", ...)
}

# ICSID: one row per tribunal seat (742 cases x 3 seats).
cases <- read("icsid_normalized.csv",
              colClasses = c(date_registered = "Date",
                             date_of_constitution_of_tribunal = "Date",
                             date_concluded = "Date"))
seats <- reshape(
  cases, direction = "long",
  varying = c("president_name", "arbitrator_1_name", "arbitrator_2_name"),
  v.names = "arbitrator", timevar = "seat",
  times = c("president", "arbitrator_1", "arbitrator_2"), idvar = "caseno"
)
icsid_tribunals <- data.frame(
  case = seats$caseno,
  seat = seats$seat,
  arbitrator = seats$arbitrator,
  registered = seats$date_registered,
  constituted = seats$date_of_constitution_of_tribunal,
  concluded = seats$date_concluded,
  is_concluded = seats$concluded == "True",
  economic_sector = seats$economic_sector,
  subject = seats$subject_of_dispute,
  respondent = seats$respondents,
  stringsAsFactors = FALSE
)
icsid_tribunals <- icsid_tribunals[order(icsid_tribunals$constituted, icsid_tribunals$case,
                                         icsid_tribunals$seat), ]
rownames(icsid_tribunals) <- NULL

# GFCC: decisions and citations, the citation dated by both decisions.
decisions <- read("bverfge_nodes_volumes-001-160.csv", colClasses = c(date = "Date"))
gfcc_decisions <- data.frame(
  decision = decisions$key, date = decisions$date, citation = decisions$bverfge,
  type = decisions$dectype, title = decisions$title, stringsAsFactors = FALSE
)
edges <- read("bverfge_edges_volumes-001-160.csv")
date_of <- stats::setNames(gfcc_decisions$date, gfcc_decisions$decision)
gfcc_citations <- data.frame(
  citing = edges$source, cited = edges$target, block = edges$chunk_id,
  date_citing = unname(date_of[edges$source]),
  date_cited = unname(date_of[edges$target]),
  stringsAsFactors = FALSE
)
gfcc_citations <- gfcc_citations[order(gfcc_citations$date_citing, gfcc_citations$citing,
                                       gfcc_citations$block, gfcc_citations$cited), ]
rownames(gfcc_citations) <- NULL

stopifnot(nrow(icsid_tribunals) == 2226L, nrow(gfcc_decisions) == 3618L,
          nrow(gfcc_citations) == 77284L, !anyNA(gfcc_citations$date_citing))
save(icsid_tribunals, file = "data/icsid_tribunals.rda", compress = "xz")
save(gfcc_decisions, file = "data/gfcc_decisions.rda", compress = "xz")
save(gfcc_citations, file = "data/gfcc_citations.rda", compress = "xz")
cat("written:", format(file.size(list.files("data", "gfcc|icsid", full.names = TRUE)) / 1024, digits = 3), "KB\n")
