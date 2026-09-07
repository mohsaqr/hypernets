test_that("contact hypergraphs are instantaneous, and cumulative on request", {
  dat <- data.frame(
    member = c("a", "b", "b", "c", "c", "d"),
    event = rep(c("e1", "e2", "e3"), each = 2),
    year = rep(1:3, each = 2),
    source = rep(c("s1", "s2", "s3"), each = 2)
  )
  thg <- temporal_hypergraph(dat, actor = "member", group = "event",
                             time = "year")
  expect_s3_class(thg, "net_temporal_hypergraph")
  expect_equal(thg$times, 1:3)
  expect_equal(thg$format, "contact")
  expect_identical(thg$time_unit, "step")
  snap_two <- hypergraph_snapshot(thg, 2)
  expect_equal(snap_two$n_hyperedges, 1)
  cumulative_two <- hypergraph_snapshot(thg, 2, mode = "cumulative")
  expect_equal(cumulative_two$n_hyperedges, 2)
  snapshots <- hypergraph_snapshots(thg)
  expect_equal(length(snapshots), 3)
  aggregate <- hypergraph_snapshot(thg, mode = "cumulative")
  expect_equal(aggregate$n_hyperedges, 3)
  memberships <- as.data.frame(thg)
  expect_equal(memberships, thg$memberships)
  # a column constant within each hyperedge is a hyperedge attribute
  expect_identical(thg$edge_data$source, c("s1", "s2", "s3"))
  sm <- summary(thg)
  expect_equal(sm$n_nodes, 4)
  expect_equal(sm$n_hyperedges, 3)
  expect_equal(sm$n_memberships, 6)
  expect_equal(sm$mean_edge_size, 2)
  expect_equal(sm$format, "contact")
})

test_that("interval snapshots use closed intervals", {
  dat <- data.frame(
    case = rep(c("A", "B"), each = 3),
    arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3"),
    start = rep(c(1, 2), each = 3), end = rep(c(3, 2), each = 3)
  )
  thg <- temporal_hypergraph(dat, actor = "arbitrator", group = "case",
                             start = "start", end = "end")
  expect_equal(thg$format, "interval")
  active_two <- hypergraph_snapshot(thg, 2)
  expect_equal(active_two$n_hyperedges, 2)
  active_three <- hypergraph_snapshot(thg, 3)
  expect_equal(active_three$n_hyperedges, 1)
  cumulative_three <- hypergraph_snapshot(thg, 3, mode = "cumulative")
  expect_equal(cumulative_three$n_hyperedges, 2)
})

test_that("an edge list gives hyperedges of size two", {
  contacts <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
                         time = 1:3, kind = c("call", "mail", "call"))
  thg <- temporal_hypergraph(contacts, from = "from", to = "to", time = "time")
  expect_equal(thg$edges, c("e1", "e2", "e3"))
  expect_equal(unname(thg$edge_data$kind), c("call", "mail", "call"))
  snap <- hypergraph_snapshot(thg, 2, mode = "cumulative")
  expect_equal(snap$n_hyperedges, 2)
  expect_equal(snap$size_distribution, c(size_2 = 2L))
  static <- group_hypergraph(contacts, from = "from", to = "to")
  expect_equal(static$n_hyperedges, 3)
  expect_equal(static$nodes, c("a", "b", "c"))
  expect_error(temporal_hypergraph(contacts, from = "from", to = "to",
                                   actor = "from", time = "time"),
               class = "hypernets_bad_input")
})

test_that("snapshot duplicate handling records multiplicity", {
  dat <- data.frame(
    member = rep(c("a", "b", "c"), 2),
    event = rep(c("e1", "e2"), each = 3), time = 1
  )
  thg <- temporal_hypergraph(dat, actor = "member", group = "event",
                             time = "time")
  multi <- hypergraph_snapshot(thg, 1, multiedges = TRUE)
  simple <- hypergraph_snapshot(thg, 1, multiedges = FALSE)
  expect_equal(multi$n_hyperedges, 2)
  expect_equal(simple$n_hyperedges, 1)
  expect_equal(simple$edge_multiplicity, 2L)
})

test_that("temporal constructors reject invalid spells and read per-row times as membership times", {
  # since 0.3.11 a time that varies within a hyperedge is a membership time:
  # the hyperedge spans the hull and each membership is a contact
  spread <- data.frame(member = c("a", "b"), event = "e1", time = c(1, 2))
  th <- temporal_hypergraph(spread, actor = "member", group = "event", time = "time")
  expect_true(th$params$membership_times)
  expect_identical(th$edge_data$start, 1)
  expect_identical(th$edge_data$end, 2)
  expect_identical(th$memberships$end, c(1, 2))
  detected <- temporal_hypergraph(spread, actor = "member", group = "event")
  expect_identical(as.data.frame(detected), as.data.frame(th))
  interval <- data.frame(member = "a", event = "e1", start = 2, end = 1)
  expect_error(temporal_hypergraph(interval, actor = "member", group = "event",
                                   start = "start", end = "end"),
               class = "hypernets_bad_input")
  expect_error(temporal_hypergraph(spread, actor = "member", group = "event",
                                   time = "time", start = "time"),
               class = "hypernets_bad_input")
})

# ---- Dynet vocabulary -------------------------------------------------------

test_that("deprecated argument names still work and warn with a class", {
  dat <- data.frame(member = c("a", "b"), event = "e1", time = 1)
  expect_warning(
    old <- temporal_hypergraph(dat, actor = "member", cooccur_by = "event", time = "time"),
    class = "hypernets_deprecated"
  )
  new <- temporal_hypergraph(dat, actor = "member", group = "event", time = "time")
  expect_identical(old, new)
  new_group <- group_hypergraph(dat, actor = "member", group = "event")
  expect_warning(
    old_member <- group_hypergraph(dat, member = "member", group = "event"),
    class = "hypernets_deprecated"
  )
  expect_warning(
    old_by <- group_hypergraph(dat, actor = "member", cooccur_by = "event"),
    class = "hypernets_deprecated"
  )
  expect_identical(old_member, new_group)
  expect_identical(old_by, new_group)
})

test_that("columns are detected from Dynet's alias table, case-insensitively", {
  log <- data.frame(Sender = c("a", "b"), Receiver = c("b", "c"),
                    Timestamp = c(1, 2), stringsAsFactors = FALSE)
  detected <- temporal_hypergraph(log)
  explicit <- temporal_hypergraph(log, from = "Sender", to = "Receiver",
                                  time = "Timestamp")
  expect_identical(detected, explicit)
  spells <- data.frame(source = c("a", "b"), target = c("b", "c"),
                       onset = c(1, 2), terminus = c(3, 4))
  interval <- temporal_hypergraph(spells)
  expect_identical(interval$format, "interval")
  expect_identical(interval$params$time, "onset")
  expect_identical(interval$params$end, "terminus")
  # an explicit name must exist as written
  expect_error(temporal_hypergraph(log, from = "sender", to = "Receiver", time = "Timestamp"),
               class = "hypernets_missing_column")
  expect_error(temporal_hypergraph(data.frame(x = 1, y = 2, time = 1)),
               class = "hypernets_bad_input")
  copresence <- data.frame(student = c("a", "b"), seminar = "s1", time = 1)
  expect_error(temporal_hypergraph(copresence, actor = "student"),
               class = "hypernets_bad_input")
  static <- group_hypergraph(spells)
  expect_identical(static$n_hyperedges, 2L)
})

test_that("calendar times become offsets from an origin in a reported unit", {
  dates <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
                      time = as.Date(c("2024-01-01", "2024-01-05", "2024-01-09")))
  as_date <- temporal_hypergraph(dates, from = "from", to = "to", time = "time")
  expect_identical(as_date$time_unit, "days")
  expect_identical(as_date$times, c(0, 4, 8))
  expect_s3_class(as_date$origin, "POSIXct")
  expect_identical(format(as_date$origin, "%Y-%m-%d"), "2024-01-01")
  # INVARIANT: character and POSIXct input give the same clock
  as_character <- temporal_hypergraph(transform(dates, time = format(time)),
                                      from = "from", to = "to", time = "time")
  as_posix <- temporal_hypergraph(transform(dates, time = as.POSIXct(time, tz = "UTC")),
                                  from = "from", to = "to", time = "time")
  expect_identical(as_character$memberships, as_date$memberships)
  expect_identical(as_posix$memberships, as_date$memberships)
  in_hours <- temporal_hypergraph(dates, from = "from", to = "to", time = "time",
                                  time_unit = "hours")
  expect_identical(in_hours$times, c(0, 96, 192))
  # a date given later is converted, and refused on a numeric clock
  snap <- hypergraph_snapshot(as_date, at = "2024-01-05")
  expect_identical(snap$params$at, 4)
  numeric <- temporal_hypergraph(transform(dates, time = c(1, 5, 9)),
                                 from = "from", to = "to", time = "time")
  expect_error(hypergraph_snapshot(numeric, at = as.Date("2024-01-05")),
               class = "hypernets_bad_input")
  expect_error(temporal_hypergraph(transform(dates, time = c("yesterday", "today", "now")),
                                   from = "from", to = "to", time = "time"),
               class = "hypernets_unparsed_time")
  expect_error(temporal_hypergraph(transform(dates, time = c(TRUE, FALSE, TRUE)),
                                   from = "from", to = "to", time = "time"),
               class = "hypernets_bad_input")
})

test_that("`time` is a contact clock and the cumulative view is a mode", {
  dat <- data.frame(member = c("a", "b", "b", "c"), event = rep(c("e1", "e2"), each = 2),
                    time = rep(c(1, 3), each = 2), stop = rep(c(2, 4), each = 2))
  expect_error(temporal_hypergraph(dat, actor = "member", group = "event",
                                   time = "time", end = "stop"),
               class = "hypernets_bad_input")
  contact <- temporal_hypergraph(dat, actor = "member", group = "event", time = "time")
  expect_identical(contact$format, "contact")
  between <- hypergraph_snapshot(contact, at = 2)
  expect_identical(between$n_hyperedges, 0L)
  cumulative <- hypergraph_snapshot(contact, at = 2, mode = "cumulative")
  expect_identical(cumulative$n_hyperedges, 1L)
  # INVARIANT: the cumulative snapshot at the end holds every hyperedge
  at_end <- hypergraph_snapshot(contact, mode = "cumulative")
  expect_identical(at_end$n_hyperedges, length(contact$edges))
  printed <- capture.output(print(contact))
  expect_true(any(grepl("contact", printed)))
  expect_true(any(grepl("cumulative", printed)))
  expect_warning(old <- hypergraph_snapshot(contact, mode = "all"),
                 class = "hypernets_deprecated")
  expect_identical(old$incidence, at_end$incidence)
})

test_that("observation bounds clip the grid without rewriting memberships", {
  seats <- data.frame(case = rep(c("A", "B"), each = 2),
                      arbitrator = c("p1", "a1", "p2", "a2"),
                      start = rep(c(1, 5), each = 2), end = rep(c(3, NA), each = 2))
  open <- temporal_hypergraph(seats, actor = "arbitrator", group = "case",
                              start = "start", end = "end")
  bounded <- temporal_hypergraph(seats, actor = "arbitrator", group = "case",
                                 start = "start", end = "end",
                                 observation_start = 2, observation_end = 8)
  expect_identical(bounded$observation, c(start = 2, end = 8))
  expect_true(bounded$params$observation_explicit)
  # INVARIANT: the stored spells are the raw ones
  expect_identical(as.data.frame(bounded), as.data.frame(open))
  expect_error(hypergraph_snapshot(bounded, at = 1), class = "hypernets_outside_observation")
  expect_error(hg_growth(bounded, start = 9), class = "hypernets_outside_observation")
  expect_error(temporal_hypergraph(seats, actor = "arbitrator", group = "case",
                                   start = "start", end = "end",
                                   observation_start = 5, observation_end = 2),
               class = "hypernets_bad_input")
  # the open-ended case is active through the end of observation
  last <- hypergraph_snapshot(bounded)
  expect_identical(last$params$at, 8)
  expect_identical(last$n_hyperedges, 1L)
  clipped <- hg_growth(bounded, start = 0, end = 10)
  expect_identical(clipped$time, c(3, 5))
  printed <- capture.output(print(bounded))
  expect_true(any(grepl("declared", printed)))
})

test_that("the measurement grid follows step, window and at", {
  contacts <- data.frame(from = c("a", "b", "c", "a"), to = c("b", "c", "a", "c"),
                         time = c(0, 1, 2, 5))
  thg <- temporal_hypergraph(contacts, from = "from", to = "to", time = "time")
  partition <- hg_growth(thg, step = 2)
  expect_identical(partition$time, c(0, 2, 4))
  # INVARIANT: a partition counts every contact exactly once
  expect_identical(sum(partition$n_edges), 4L)
  rolling <- hg_growth(thg, step = 1, window = 3)
  expect_identical(rolling$n_edges[1L], 3L)
  points <- hg_growth(thg, step = 1, window = 0)
  expect_identical(points$n_edges, c(1L, 1L, 1L, 0L, 0L, 1L))
  whole <- hypergraph_snapshot(thg, window = "all")
  expect_identical(whole$n_hyperedges, 4L)
  snaps <- hypergraph_snapshots(thg, step = 2)
  expect_identical(names(snaps), c("0", "2", "4"))
  expect_identical(attr(snaps, "window"), 2)
  edges <- hg_edges(thg, step = 2)
  expect_identical(unique(edges$time), c(0, 2, 4))
  instants <- hypergraph_snapshots(thg, at = c(5, 1))
  expect_identical(names(instants), c("1", "5"))
  expect_error(hg_growth(thg, step = 2, window = "all"), class = "hypernets_bad_input")
  expect_error(hg_growth(thg, at = 1, step = 1), class = "hypernets_bad_input")
  expect_error(hg_growth(thg, step = 0), class = "hypernets_bad_input")
  expect_error(hg_growth(thg, window = -1), class = "hypernets_bad_input")
  expect_error(hg_growth(thg, start = 3, end = 1), class = "hypernets_bad_input")
  expect_error(hypergraph_snapshot(thg, at = c(1, 2)), class = "hypernets_bad_input")
})

test_that("INVARIANT: hypernets and Dynet agree on which actor pairs are ever co-present", {
  skip_if_not_installed("Dynet")
  log <- data.frame(
    student = c("a", "b", "c", "a", "b", "d", "c", "d"),
    seminar = c("s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3"),
    start = as.Date(c(rep("2024-01-01", 3), rep("2024-02-01", 3), rep("2024-03-01", 2))),
    end = as.Date(c(rep("2024-01-10", 3), rep("2024-02-05", 3), rep("2024-03-02", 2))),
    stringsAsFactors = FALSE
  )
  dn <- Dynet::dynet(log, actor = "student", group = "seminar")
  dynet_spells <- as.data.frame(dn)
  dynet_pairs <- unique(paste(pmin(dynet_spells$from, dynet_spells$to),
                              pmax(dynet_spells$from, dynet_spells$to)))
  thg <- temporal_hypergraph(log, actor = "student", group = "seminar",
                             start = "start", end = "end")
  aggregate <- hypergraph_snapshot(thg, mode = "cumulative")
  projection <- hg_project(aggregate)
  hypernets_pairs <- unique(paste(pmin(projection$from, projection$to),
                               pmax(projection$from, projection$to)))
  expect_setequal(hypernets_pairs, dynet_pairs)
  # and on the clock: the same unit and origin
  expect_identical(thg$time_unit, dn$meta$time_unit)
  expect_identical(as.numeric(thg$origin), as.numeric(dn$meta$origin))
  expect_identical(sort(unique(thg$memberships$start)), sort(unique(dynet_spells$start)))
})
