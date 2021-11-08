#' @include utils.R
NULL

#' Pre-Defined Time-to-Event Source Objects
#'
#' @export
#'
#' @rdname tte_source_objects
#'
death_event <- event_source(
  dataset_name = "adsl",
  filter = DTHFL == "Y",
  date = DTHDT,
  set_values_to = vars(
    EVNTDESC = "DEATH",
    SRCDOM = "ADSL",
    SRCVAR = "DTHDT"
  )
)

#' @rdname tte_source_objects
#' @export
lastalive_censor <- censor_source(
  dataset_name = "adsl",
  date = LSTALVDT,
  set_values_to = vars(
    EVNTDESC = "ALIVE",
    SRCDOM = "ADSL",
    SRCVAR = "LSTALVDT"
  )
)

#' @rdname tte_source_objects
#' @export
ae_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y",
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "ADVERSE EVENT",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)

#' @rdname tte_source_objects
#' @export
ae_ser_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y" & AESER == "Y" ,
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "SERIOUS ADVERSE EVENT",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)

#' @rdname tte_source_objects
#' @export
ae_gr1_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y" & ATOXGR == "1",
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "GRADE 1 ADVERSE EVENT",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)

#' @rdname tte_source_objects
#' @export
ae_gr2_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y" & ATOXGR == "2",
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "GRADE 2 ADVERSE EVENT",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)

#' @rdname tte_source_objects
#' @export
ae_gr3_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y" & ATOXGR == "3",
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "GRADE 3 ADVERSE EVENT",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)

#' @rdname tte_source_objects
#' @export
ae_gr4_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y" & ATOXGR == "4",
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "GRADE 4 ADVERSE EVENT",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)

#' @rdname tte_source_objects
#' @export
ae_gr5_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y" & ATOXGR == "5",
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "GRADE 5 ADVERSE EVENT",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)

#' @rdname tte_source_objects
#' @export
ae_gr35_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y" & ATOXGR %in% c("3", "4", "5"),
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "GRADE 5 ADVERSE EVENT",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)

#' @rdname tte_source_objects
#' @export
ae_sev_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y" & AESEV == "SEVERE",
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "SEVERE ADVERSE EVENT",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)

#' @rdname tte_source_objects
#' @export
ae_wd_event <- event_source(
  dataset_name = "adae",
  filter = TRTEMFL == "Y" & AEACN == "DRUG WITHDRAWN",
  date = ASTDT,
  set_values_to = vars(
    EVNTDESC = "ADVERSE EVENT LEADING TO DRUG WITHDRAWAL",
    SRCDOM = "ADAE",
    SRCVAR = "ASTDT",
    SRCSEQ = AESEQ
  )
)
