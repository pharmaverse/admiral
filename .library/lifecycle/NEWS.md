
# lifecycle 0.1.0

* Deprecated functions under the control of the developer now warn
  repeatedly in unit tests.

* Deprecation warnings now record a backtrace. Call
  `lifecycle::last_warnings()` and `lifecycle::last_warning()` to
  print the warnings that occurred during the last command, along with
  their backtraces.

* The naming scheme of signaller functions has been simplified:

  - `signal_soft_deprecated()` is now `deprecate_soft()`.
  - `warn_deprecated()` is now `deprecate_warn()`.
  - `stop_defunct()` is now `deprecate_stop()`.

* The signaller functions now take a version and two descriptors for
  the deprecated feature and its replacement (the latter is
  optional). The deprecation message is built from these
  components. You can pass a `details` argument to append additional
  information to the generated deprecation message.

* Helpers from rlang's `compat-lifecycle.R` drop-in file are now
  exported in this package.
