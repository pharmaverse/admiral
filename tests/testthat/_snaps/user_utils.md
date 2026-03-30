# print.source Test 13: `source` objects are printed as intended

    Code
      print(ttae)
    Output
      <event_source> object
      dataset_name: "ae"
      filter: NULL
      date: AESTDTC
      censor: 0
      set_values_to:
        EVENTDESC: "AE"
        SRCDOM: "AE"
        SRCVAR: "AESTDTC"
        SRCSEQ: AESEQ
      order: NULL
      consider_end_dates: TRUE

# print_named_list Test 18: named list with unamed list

    Code
      print_named_list(list(list_item = list("Hello World!", expr(universe), list(42)),
      another_one = ymd("2020-02-02")))
    Output
      list_item:
        "Hello World!"
        universe
        42
      another_one: 2020-02-02

