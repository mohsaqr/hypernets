# Build the bundled `ring_sequences` and `ring_communities` datasets: a
# simulated process with planted memory modules, used by the
# hon_communities() documentation and tests.
# Run from the package root: Rscript data-raw/ring_sequences.R
#
# Four groups of four actions sit on a ring; neighbouring groups share one
# action (s1..s4), the other two actions of a group (p<k>_1, p<k>_2) are its
# own. A walk starts in a random group at a random action of it; each step
# moves to another action of the current group, chosen at random; on
# reaching a shared action the walk switches to the neighbouring group with
# probability 0.05. 200 walks of 40 steps, seed 1.

set.seed(1L)
n_sequences <- 200L
walk_length <- 40L
n_groups <- 4L
group_size <- 4L
switch_probability <- 0.05

shared <- sprintf("s%d", seq_len(n_groups))
groups <- lapply(seq_len(n_groups), \(k) {
  c(shared[k], sprintf("p%d_%d", k, seq_len(group_size - 2L)),
    shared[k %% n_groups + 1L])
})

# A walk is sequential (each step starts where the last ended, in the group
# the last one left it in), so the steps are carried by Reduce(); the
# random draws are made in walk order, one walk after another.
one_walk <- function(i) {
  start_group <- sample.int(n_groups, 1L)
  start <- list(action = sample(groups[[start_group]], 1L), group = start_group)
  steps <- Reduce(\(state, step) {
    action <- sample(setdiff(groups[[state$group]], state$action), 1L)
    holding <- which(vapply(groups, \(g) action %in% g, logical(1L)))
    other <- setdiff(holding, state$group)
    switch_now <- length(other) > 0L && stats::runif(1L) < switch_probability
    list(action = action, group = if (switch_now) other else state$group)
  }, seq.int(2L, walk_length), start, accumulate = TRUE)
  vapply(steps, \(state) state$action, character(1L))
}
ring_sequences <- lapply(seq_len(n_sequences), one_walk)

# The planted community of every state a second-order network can hold: a
# state "u -> v" belongs to the one group holding both u and v, and the
# first-order state of an action that is not shared to that action's group.
# First-order states of the shared actions have no single group and are
# left out.
pairs <- do.call(rbind, lapply(seq_len(n_groups), \(k) {
  both <- expand.grid(from = groups[[k]], to = groups[[k]],
                      stringsAsFactors = FALSE)
  both <- subset(both, from != to)
  data.frame(node = paste(both$from, both$to, sep = " -> "), community = k)
}))
private <- data.frame(
  node = unlist(lapply(groups, \(g) setdiff(g, shared))),
  community = rep(seq_len(n_groups), each = group_size - 2L)
)
ring_communities <- rbind(private, pairs)
ring_communities <- ring_communities[order(ring_communities$community,
                                           ring_communities$node), ]
rownames(ring_communities) <- NULL
ring_communities$community <- as.integer(ring_communities$community)

stopifnot(
  "one sequence per walk" = length(ring_sequences) == n_sequences,
  "every walk is 40 steps" = all(lengths(ring_sequences) == walk_length),
  "no node has two communities" = !anyDuplicated(ring_communities$node),
  "8 private + 4 x 12 ordered pairs" = nrow(ring_communities) == 56L
)

save(ring_sequences, file = file.path("data", "ring_sequences.rda"),
     compress = "xz")
save(ring_communities, file = file.path("data", "ring_communities.rda"),
     compress = "xz")
