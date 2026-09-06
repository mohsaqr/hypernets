# Bundled example data, shared with Nestimate (same .rda files, same
# provenance); bundled here so higher-order examples, inference, and
# tutorials run on real coded sequences.

#' Human-AI Vibe Coding Interaction Data (Long Format)
#'
#' Coded turns from 429 human-AI pair programming sessions across 34
#' projects, in long format: `human_long` holds the human turns (10,796
#' rows), `ai_long` the AI turns (8,551 rows). Each session's ordered
#' codes form one categorical sequence, which makes the pair a natural
#' two-cohort input for the higher-order verbs -- e.g.
#' `bootstrap_hon(human_long, action = "code", actor = "session_id",
#' time = "timestamp")` or `compare_hon(human_long, ai_long, ...)`.
#'
#' The same data feed all three structure families: the ordered codes are
#' sequences for the memory family, a session is a natural hyperedge over
#' the codes that co-occur in it ([group_hypergraph()]), and a fitted
#' memory network becomes a pathway complex for the simplicial family
#' ([build_simplicial()] with `type = "pathway"`).
#'
#' @format Data frames in long format with 9 columns:
#' \describe{
#'   \item{message_id}{Integer. Turn index.}
#'   \item{project}{Character. Project identifier (Project_1 .. Project_34).}
#'   \item{session_id}{Character. Unique session hash.}
#'   \item{timestamp}{Integer. Unix timestamp for ordering.}
#'   \item{session_date}{Character. Date of the session (YYYY-MM-DD).}
#'   \item{code}{Character. Interaction code.}
#'   \item{cluster}{Character. High-level cluster: Directive, Evaluative,
#'     or Metacognitive (human codes); AI turns carry their own scheme.}
#'   \item{code_order}{Integer. Order of the code within the session.}
#'   \item{order_in_session}{Integer. Absolute turn order within the session.}
#' }
#'
#' @source Saqr, M. (2026). Human-AI vibe coding interaction study.
#'   \url{https://saqr.me/blog/2026/human-ai-interaction-cograph/}
#'
#' @examples
#' bs <- bootstrap_hon(human_long, action = "code", actor = "session_id",
#'                     time = "timestamp", n_boot = 20, max_order = 2,
#'                     seed = 1)
#' rules <- as.data.frame(bs, order_min = 2)
#' head(rules)
#'
#' @name long-data
#' @aliases ai_long
NULL

#' @rdname long-data
"human_long"

#' @rdname long-data
"ai_long"

#' ICSID arbitration tribunals, one row per seat
#'
#' The 742 cases registered at the International Centre for Settlement of
#' Investment Disputes between 1974 and June 2023 for which the three
#' tribunal members and the dates are known, as released by Coupette,
#' Hartung and Katz (2024): one row per seat, so a case has three rows.
#' Arbitrators sharing a case are the co-occurrence data of a temporal
#' hypergraph whose hyperedges are tribunals active from constitution to
#' conclusion; see `vignette("legal-hypergraphs")`.
#'
#' @format A data frame with 2,226 rows and 10 columns:
#' \describe{
#'   \item{case}{ICSID case number, the hyperedge.}
#'   \item{seat}{`"president"`, `"arbitrator_1"` or `"arbitrator_2"`.}
#'   \item{arbitrator}{Name of the arbitrator holding the seat, the node.}
#'   \item{registered}{Date the case was registered.}
#'   \item{constituted}{Date the tribunal was constituted.}
#'   \item{concluded}{Date the case concluded; `NA` for a pending case, whose
#'     tribunal stays active through the end of observation (2023-06-15).}
#'   \item{is_concluded}{Whether the archive marks the case as concluded.}
#'   \item{economic_sector}{Economic sector of the dispute.}
#'   \item{subject}{Subject of the dispute.}
#'   \item{respondent}{Respondent state.}
#' }
#' @source Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs
#'   (ICSID dataset). Zenodo. \doi{10.5281/zenodo.8081513}, via the
#'   replication archive \doi{10.5281/zenodo.8081507}. Licence CC BY-NC 4.0.
#'   Built by `data-raw/legal_hypergraphs.R`.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' head(icsid_tribunals)
#' tribunals <- temporal_hypergraph(icsid_tribunals, actor = "arbitrator",
#'                                  group = "case", start = "constituted",
#'                                  end = "concluded")
#' summary(tribunals)
"icsid_tribunals"

#' Decisions of the German Federal Constitutional Court
#'
#' The 3,618 decisions in volumes 1 to 160 of the official collection of the
#' German Federal Constitutional Court (BVerfGE), 1951 to 2022, as released
#' by Coupette, Hartung and Katz (2024). The node table of the citation-block
#' hypergraph in [gfcc_citations]: every decision is a node from its own
#' date, whether or not it is ever cited.
#'
#' @format A data frame with 3,618 rows and 5 columns:
#' \describe{
#'   \item{decision}{Decision key, `volume-page`.}
#'   \item{date}{Date of the decision.}
#'   \item{citation}{Official citation, e.g. `"BVerfGE 1, 1-3"`.}
#'   \item{type}{Decision type (`"Urteil"` or `"Beschluss"`).}
#'   \item{title}{Title of the decision, in German.}
#' }
#' @source Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs
#'   (GFCC dataset). Zenodo. \doi{10.5281/zenodo.8081511}, via the
#'   replication archive \doi{10.5281/zenodo.8081507}. Licence CC BY-NC 4.0.
#'   Built by `data-raw/legal_hypergraphs.R`.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' head(gfcc_decisions)
"gfcc_decisions"

#' Citation blocks of the German Federal Constitutional Court
#'
#' Every citation between the decisions in [gfcc_decisions], grouped into
#' citation blocks: an uninterrupted run of cited decisions inside one
#' decision, the hyperedge of Coupette, Hartung and Katz (2024). One row per
#' citation, dated by both decisions, so backward citations are the rows
#' with `date_citing > date_cited`.
#'
#' @format A data frame with 77,284 rows and 5 columns:
#' \describe{
#'   \item{citing}{Key of the citing decision, the source of the block.}
#'   \item{cited}{Key of the cited decision, the node.}
#'   \item{block}{Citation block identifier, `citing:index`, the hyperedge.}
#'   \item{date_citing}{Date of the citing decision, the time of the block.}
#'   \item{date_cited}{Date of the cited decision.}
#' }
#' @source Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs
#'   (GFCC dataset). Zenodo. \doi{10.5281/zenodo.8081511}, via the
#'   replication archive \doi{10.5281/zenodo.8081507}. Licence CC BY-NC 4.0.
#'   Built by `data-raw/legal_hypergraphs.R`.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' head(gfcc_citations)
#' blocks <- temporal_hypergraph(gfcc_citations, actor = "cited",
#'                               group = "block", time = "date_citing",
#'                               nodes = gfcc_decisions, sparse = TRUE)
#' summary(blocks)
"gfcc_citations"
