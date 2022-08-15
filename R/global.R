# R CMD check will generate NOTE(s) with message 'no visible binding for
# global variable xyz'. Adding this 'xyz' variable to this file will remove
# this note. See `?globalVariables` for more information.
globalVariables(c(
  "_unit",
  "auto",
  "PARAMCD"
))
