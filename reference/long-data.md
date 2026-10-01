# Human-AI Vibe Coding Interaction Data (Long Format)

Coded turns from 429 human-AI pair programming sessions across 34
projects, in long format: `human_long` holds the human turns (10,796
rows), `ai_long` the AI turns (8,551 rows). Each session's ordered codes
form one categorical sequence, which makes the pair a natural two-cohort
input for the higher-order verbs – e.g.
`hg_bootstrap(human_long, action = "code", actor = "session_id", time = "timestamp")`
or `hg_compare(human_long, ai_long, ...)`.

## Usage

``` r
human_long

ai_long
```

## Format

Data frames in long format with 9 columns:

- message_id:

  Integer. Turn index.

- project:

  Character. Project identifier (Project_1 .. Project_34).

- session_id:

  Character. Unique session hash.

- timestamp:

  Integer. Unix timestamp for ordering.

- session_date:

  Character. Date of the session (YYYY-MM-DD).

- code:

  Character. Interaction code.

- cluster:

  Character. High-level cluster: Directive, Evaluative, or Metacognitive
  (human codes); AI turns carry their own scheme.

- code_order:

  Integer. Order of the code within the session.

- order_in_session:

  Integer. Absolute turn order within the session.

An object of class `data.frame` with 10796 rows and 9 columns.

An object of class `data.frame` with 8551 rows and 9 columns.

## Source

Saqr, M. (2026). Human-AI vibe coding interaction study.
<https://saqr.me/blog/2026/human-ai-interaction-cograph/>

## Details

The same data feed all three structure families: the ordered codes are
sequences for the memory family, a session is a natural hyperedge over
the codes that co-occur in it
([`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md)),
and a fitted memory network becomes a pathway complex for the simplicial
family
([`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md)
with `type = "pathway"`).

## Examples

``` r
bs <- hg_bootstrap(human_long, action = "code", actor = "session_id",
                    time = "timestamp", n_boot = 20, max_order = 2,
                    seed = 1)
#> Error in Nestimate::prepare(data[columns], actor = actor, action = action,     time = time, session = session, time_threshold = time_threshold,     timezone = timezone): unused argument (timezone = timezone)
rules <- hg_get(bs, order_min = 2)
#> Error: object 'bs' not found
head(rules)
#> Error: object 'rules' not found
```
