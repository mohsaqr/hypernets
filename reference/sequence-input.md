# Sequence input shared by the sequence verbs

Sequence input shared by the sequence verbs

## Sequence input

Every verb that reads sequences reads them the same way:

- Long event table:

  a data.frame with one event per row. `action` names the state column,
  `actor` the column (or columns) that identify whose events they are,
  `session` the column (or columns) that split an actor's events into
  sessions, and `time` the column that orders them: timestamps in any
  common date-time format, Unix times, or positions; row order when
  `NULL`. A column named `action`, `time`, `session` or `session_id` (in
  any case) is used when its argument is `NULL`; `session = FALSE`
  switches the session detection off. When `time` is given, a gap of
  more than `time_threshold` seconds between consecutive events starts a
  new sequence (default `900`, fifteen minutes);
  `time_threshold = FALSE` keeps each actor or session in one sequence
  however long the gaps. `timezone` is the time zone of timestamps that
  carry none (default `"UTC"`). Events with a missing actor or session
  raise `hypernets_bad_input`; a missing action stays in its sequence as
  a gap. Without `actor` all events form one sequence, announced by the
  message `hypernets_single_sequence`.

- Wide data.frame or character matrix:

  one sequence per row; trailing `NA`s end a sequence.

- List:

  one character vector per sequence.

- Model object:

  a `netobject`, `netobject_group`, `tna` or `cograph_network` that
  carries its sequence data.

These are the conventions of the tna family of packages, and the same
call builds the same sequences there.

A data.frame with columns named like an event table (`code`, `state`,
`user`, `timestamp`, ...) that is passed without `action =` and has no
`action` column raises `hypernets_long_format` (a
`hypernets_bad_input`): read as wide, its actor ids and times would
silently become states. `actor`, `time` or `session` without an action
column raises `hypernets_bad_input`.
