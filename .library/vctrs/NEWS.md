
# vctrs 0.2.3

* The main feature of this release is considerable performance
  improvements with factors and dates.

* `vec_c()` now falls back to `base::c()` if the vector doesn't
  implement `vec_ptype2()` but implements `c()`. This should improve
  the compatibility of vctrs-based functions with foreign classes
  (#801).

* `new_data_frame()` is now faster.

* New `vec_is_list()` for detecting if a vector is a list in the vctrs sense.
  For instance, objects of class `lm` are not lists. In general, classes need
  to explicitly inherit from `"list"` to be considered as lists by vctrs.

* Unspecified vectors of `NA` can now be assigned into a list (#819).

  ```
  x <- list(1, 2)
  vec_slice(x, 1) <- NA
  x
  #> [[1]]
  #> NULL
  #>
  #> [[2]]
  #> 2
  ```

* `vec_ptype()` now errors on scalar inputs (#807).

* `vec_ptype_finalise()` is now recursive over all data frame types, ensuring
  that unspecified columns are correctly finalised to logical (#800).

* `vec_ptype()` now correctly handles unspecified columns in data frames, and
  will always return an unspecified column type (#800).

* `vec_slice()` and `vec_chop()` now work correctly with `bit64::integer64()`
  objects when an `NA` subscript is supplied. By extension, this means that
  `vec_init()` now works with these objects as well (#813).

* `vec_rbind()` now binds row names. When named inputs are supplied
  and `names_to` is `NULL`, the names define row names. If `names_to`
  is supplied, they are assigned in the column name as before.

* `vec_cbind()` now binds row names if they are congruent across
  inputs. If the row names are not identical that's an error.

* The `c()` method for `vctrs_vctr` now throws an error when
  `recursive` or `use.names` is supplied (#791).


# vctrs 0.2.2

* New `vec_as_subscript()` function to cast inputs to the base type
  of a subscript (logical, numeric, or character). `vec_as_index()`
  has been renamed to `vec_as_location()`. Use `num_as_location()` if
  you need more options to control how numeric subscripts are
  converted to a vector of locations.

* New `vec_as_subscript2()`, `vec_as_location2()`, and
  `num_as_location2()` variants for validating scalar subscripts and
  locations (e.g. for indexing with `[[`).

* `vec_as_location()` now preserves names of its inputs if possible.

* `vec_ptype2()` methods for base classes now prevent
  inheritance. This makes sense because the subtyping graph created by
  `vec_ptype2()` methods is generally not the same as the inheritance
  relationships defined by S3 classes. For instance, subclasses are
  often a richer type than their superclasses, and should often be
  declared as supertypes (e.g. `vec_ptype2()` should return the
  subclass).

  We introduced this breaking change in a patch release because
  `new_vctr()` now adds the base type to the class vector by default,
  which caused `vec_ptype2()` to dispatch erroneously to the methods
  for base types. We'll finish switching to this approach in vctrs
  0.3.0 for the rest of the base S3 classes (dates, data frames, ...).

* `vec_equal_na()` now works with complex vectors.

* `vctrs_vctr` class gains an `as.POSIXlt()` method (#717).

* `vec_is()` now ignores names and row names (#707).

* `vec_slice()` now support Altvec vectors (@jimhester, #696).

* `vec_proxy_equal()` is now applied recursively across the columns of
  data frames (#641).

* `vec_split()` no longer returns the `val` column as a `list_of`. It is now
  returned as a bare list (#660).

* Complex numbers are now coercible with integer and double (#564).

* zeallot has been moved from Imports to Suggests, meaning that `%<-%` is no
  longer re-exported from vctrs.

* `vec_equal()` no longer propagates missing values when comparing list
  elements. This means that `vec_equal(list(NULL), list(NULL))` will continue to
  return `NA` because `NULL` is the missing element for a list, but now
  `vec_equal(list(NA), list(NA))` returns `TRUE` because the `NA` values are
  compared directly without checking for missingness.

* Lists of expressions are now supported in `vec_equal()` and functions that
  compare elements, such as `vec_unique()` and `vec_match()`. This ensures that
  they work with the result of modeling functions like `glm()` and `mgcv::gam()`
  which store "family" objects containing expressions (#643).

* `new_vctr()` gains an experimental `inherit_base_type` argument
  which determines whether or not the class of the underlying type
  will be included in the class.

* `list_of()` now inherits explicitly from "list" (#593).

* `vec_ptype()` has relaxed default behaviour for base types; now if two
  vectors both inherit from (e.g.) "character", the common type is also
  "character" (#497).

* `vec_equal()` now correctly treats `NULL` as the missing value element for
  lists (#653).

* `vec_cast()` now casts data frames to lists rowwise, i.e. to a list of
  data frames of size 1. This preserves the invariant of
  `vec_size(vec_cast(x, to)) == vec_size(x)` (#639).

* Positive and negative 0 are now considered equivalent by all functions that
  check for equality or uniqueness (#637).

* New experimental functions `vec_group_rle()` for returning run
  length encoded groups; `vec_group_id()` for constructing group
  identifiers from a vector; `vec_group_loc()` for computing the
  locations of unique groups in a vector (#514).

* New `vec_chop()` for repeatedly slicing a vector. It efficiently captures
  the pattern of `map(indices, vec_slice, x = x)`.

* Support for multiple character encodings has been added to functions that
  compare elements within a single vector, such as `vec_unique()`, and across
  multiple vectors, such as `vec_match()`. When multiple encodings are
  encountered, a translation to UTF-8 is performed before any comparisons are
  made (#600, #553).

* Equality and ordering methods are now implemented for raw and
  complex vectors (@romainfrancois).


# vctrs 0.2.1

Maintenance release for CRAN checks.


# vctrs 0.2.0

With the 0.2.0 release, many vctrs functions have been rewritten with
native C code to improve performance. Functions like `vec_c()` and
`vec_rbind()` should now be fast enough to be used in packages. This
is an ongoing effort, for instance the handling of factors and dates
has not been rewritten yet. These classes still slow down vctrs
primitives.

The API in 0.2.0 has been updated, please see a list of breaking
changes below. vctrs has now graduated from experimental to a maturing
package (see the [lifecycle of tidyverse packages](https://www.tidyverse.org/lifecycle/)).
Please note that API changes are still planned for future releases,
for instance `vec_ptype2()` and `vec_cast()` might need to return a
sentinel instead of failing with an error when there is no common type
or possible cast.


## Breaking changes

* Lossy casts now throw errors of type `vctrs_error_cast_lossy`.
  Previously these were warnings. You can suppress these errors
  selectively with `allow_lossy_cast()` to get the partial cast
  results. To implement your own lossy cast operation, call the new
  exported function `maybe_lossy_cast()`.

* `vec_c()` now fails when an input is supplied with a name but has
  internal names or is length > 1:

  ```
  vec_c(foo = c(a = 1))
  #> Error: Can't merge the outer name `foo` with a named vector.
  #> Please supply a `.name_spec` specification.

  vec_c(foo = 1:3)
  #> Error: Can't merge the outer name `foo` with a vector of length > 1.
  #> Please supply a `.name_spec` specification.
  ```

  You can supply a name specification that describes how to combine
  the external name of the input with its internal names or positions:

  ```
  # Name spec as glue string:
  vec_c(foo = c(a = 1), .name_spec = "{outer}_{inner}")

  # Name spec as a function:
  vec_c(foo = c(a = 1), .name_spec = function(outer, inner) paste(outer, inner, sep = "_"))
  vec_c(foo = c(a = 1), .name_spec = ~ paste(.x, .y, sep = "_"))
  ```

* `vec_empty()` has been renamed to `vec_is_empty()`.

* `vec_dim()` and `vec_dims()` are no longer exported.

* `vec_na()` has been renamed to `vec_init()`, as the primary use case
  is to initialize an output container.

* `vec_slice<-` is now type stable (#140). It always returns the same
  type as the LHS. If needed, the RHS is cast to the correct type, but
  only if both inputs are coercible. See examples in `?vec_slice`.

* We have renamed the `type` particle to `ptype`:

  - `vec_type()` => `vec_ptype()`
  - `vec_type2()` => `vec_ptype2()`
  - `vec_type_common()` => `vec_ptype_common()`

  Consequently, `vec_ptype()` was renamed to `vec_ptype_show()`.


## New features

* New `vec_proxy()` generic. This is the main customisation point in
  vctrs along with `vec_restore()`. You should only implement it when
  your type is designed around a non-vector class (atomic vectors,
  bare lists, data frames). In this case, `vec_proxy()` should return
  such a vector class. The vctrs operations will be applied on the
  proxy and `vec_restore()` is called to restore the original
  representation of your type.

  The most common case where you need to implement `vec_proxy()` is
  for S3 lists. In vctrs, S3 lists are treated as scalars by
  default. This way we don't treat objects like model fits as
  vectors. To prevent vctrs from treating your S3 list as a scalar,
  unclass it from the `vec_proxy()` method. For instance here is the
  definition for `list_of`:

  ```
  #' @export
  vec_proxy.vctrs_list_of <- function(x) {
    unclass(x)
  }
  ```

  If you inherit from `vctrs_vctr` or `vctrs_rcrd` you don't need to
  implement `vec_proxy()`.

* `vec_c()`, `vec_rbind()`, and `vec_cbind()` gain a `.name_repair`
  argument (#227, #229).

* `vec_c()`, `vec_rbind()`, `vec_cbind()`, and all functions relying
  on `vec_ptype_common()` now have more informative error messages
  when some of the inputs have nested data frames that are not
  convergent:

  ```
  df1 <- tibble(foo = tibble(bar = tibble(x = 1:3, y = letters[1:3])))
  df2 <- tibble(foo = tibble(bar = tibble(x = 1:3, y = 4:6)))

  vec_rbind(df1, df2)
  #> Error: No common type for `..1$foo$bar$y` <character> and `..2$foo$bar$y` <integer>.
  ```

* `vec_cbind()` now turns named data frames to packed columns.

  ```r
  data <- tibble::tibble(x = 1:3, y = letters[1:3])
  data <- vec_cbind(data, packed = data)
  data
  # A tibble: 3 x 3
        x y     packed$x $y
    <int> <chr>    <int> <chr>
  1     1 a            1 a
  2     2 b            2 b
  3     3 c            3 c
  ```

  Packed data frames are nested in a single column. This makes it
  possible to access it through a single name:

  ```r
  data$packed
  # A tibble: 3 x 2
        x y
    <int> <chr>
  1     1 a
  2     2 b
  3     3 c
  ```

  We are planning to use this syntax more widely in the tidyverse.

* New `vec_is()` function to check whether a vector conforms to a
  prototype and/or a size. Unlike `vec_assert()`, it doesn't throw
  errors but returns `TRUE` or `FALSE` (#79).

  Called without a specific type or size, `vec_assert()` tests whether
  an object is a data vector or a scalar. S3 lists are treated as
  scalars by default. Implement a `vec_is_vector()` for your class to
  override this property (or derive from `vctrs_vctr`).

* New `vec_order()` and `vec_sort()` for ordering and sorting
  generalised vectors.

* New `.names_to` parameter for `vec_rbind()`. If supplied, this
  should be the name of a column where the names of the inputs are
  copied. This is similar to the `.id` parameter of
  `dplyr::bind_rows()`.

* New `vec_seq_along()` and `vec_init_along()` create useful sequences (#189).

* `vec_slice()` now preserves character row names, if present.

* New `vec_split(x, by)` is a generalisation of `split()` that can divide
  a vector into groups formed by the unique values of another vector. Returns
  a two-column data frame containing unique values of `by` aligned with
  matching `x` values (#196).


## Other features and bug fixes

* Using classed errors of class `"vctrs_error_assert"` for failed
  assertions, and of class `"vctrs_error_incompatible"` (with
  subclasses `_type`, `_cast` and `_op`) for errors on incompatible
  types (#184).

* Character indexing is now only supported for named objects, an error
  is raised for unnamed objects (#171).

* Predicate generics now consistently return logical vectors when
  passed a `vctrs_vctr` class. They used to restore the output to
  their input type (#251).

* `list_of()` now has an `as.character()` method. It uses
  `vec_ptype_abbr()` to collapse complex objects into their type
  representation (tidyverse/tidyr#654).

* New `stop_incompatible_size()` to signal a failure due to mismatched sizes.

* New `validate_list_of()` (#193).

* `vec_arith()` is consistent with base R when combining `difftime`
  and `date`, with a warning if casts are lossy (#192).

* `vec_c()` and `vec_rbind()` now handle data.frame columns properly
  (@yutannihilation, #182).

* `vec_cast(x, data.frame())` preserves the number of rows in `x`.

* `vec_equal()` now handles missing values symmetrically (#204).

* `vec_equal_na()` now returns `TRUE` for data frames and records when
  every component is missing, not when _any_ component is missing
  (#201).

* `vec_init()` checks input is a vector.

* `vec_proxy_compare()` gains an experimental `relax` argument, which
  allows data frames to be orderable even if all their columns are not
  (#210).

* `vec_size()` now works with positive short row names. This fixes
  issues with data frames created with jsonlite (#220).

* `vec_slice<-` now has a `vec_assign()` alias. Use `vec_assign()`
  when you don't want to modify the original input.

* `vec_slice()` now calls `vec_restore()` automatically. Unlike the
  default `[` method from base R, attributes are preserved by default.

* `vec_slice()` can correct slice 0-row data frames (#179).

* New `vec_repeat()` for repeating each element of a vector the same number
  of times.

* `vec_type2(x, data.frame())` ensures that the returned object has
  names that are a length-0 character vector.
